%Version 1: Only one file at a time

%Returns the indices of the start and stop of a sharp wave oscillation along
%with the relative index of the maximum oscillation power. It can also
%return the 
%
%
%150 300 15 500 50
function [oscIdxVec, varargout] = getOscillationIndex(interestRawVec, baselineRawVec, highPass, lowPass, minLength, maxLength, minGap, varargin)
    
    oscIdxVec = [];
    bigPowEnvelope = [];
    bigPow = 0;
    bigPowChan = 0;

    %default variables
    varStrings = ["fs" "verbose" "develop" "swrMaxLength" "swrMinLength" "swrMinGap" "correctSWR" "strongThresh" "weakThresh" "chanExamine" "passType"];
    fs = 24414;
    correctIdx = 1; %If 1 the program will go through and combine any oscillation events that have an intraoscillation time of less than or equal to the indicated gap length. With anything remaining, it'll then get rid of any SWR that are shorter than the minimum length or longer than the maximum length.
    verbose = 1;
    develop = 0;
    peakPow = 0; %-1 only looks at the negative deflecting waves for peak power, 0 looks at both, and 1 looks at the postive deflecting waves.
    strongThresh = 7; %STD
    weakThresh = 0.5; %STD
    chanExamine = []; %Which channels should be looked at. An empty matrix means all of them.
    passType = 'bandpass';
    for i = 1:2:length(varargin)
        if ~ismember(varargin{i}, varStrings) 
            fprintf('\n\nInput does not match allowable options.\nYou entered %s\nPlease try again.\n', string(varargin{i}))
            return;
        end
        eval([varargin{i} '=varargin{i + 1};']);
    end
    if develop; verbose = 1; end
    maxIdx = floor(fs / 1000); %Starting a millisecond into recording to get rid of edge cases.
    
    if passType == 'bandpass'
        [b, a] = butter(4, [highPass lowPass] ./ (fs / 2), passType);
    elseif passType == 'lowpass'
        [b, a] = butter(4, [lowPass] / (fs / 2), passType);
    else
        [b, a] = butter(4, [highPass] / (fs / 2), passType);
    end
    
    if ~isempty(baselineRawVec)
        if ~isnumeric(baselineRawVec)
            if verbose; fprintf('\nReading in baseline raw data'); end
            baselineRawVec = readmda(char(baselineRawVec));
        end
        baselineRawVec = double(baselineRawVec);

        if isempty(chanExamine); chanExamine = 1:size(baselineRawVec, 1); end

        fprintf('\n\nExamining baseline epoch to determine oscillation thresholds.')
        fprintf('\nFinding channel with largest oscillation power');
        for chan = chanExamine
            fprintf('\nAnalysing channel: %d', chan)

            rawFilt = filtfilt(b, a, baselineRawVec(chan, :));
            rawPow = rawFilt.^2;
            clear rawFilt;
            tempTotPow = sum(rawPow);
            if tempTotPow > bigPow
                bigPow = tempTotPow;
                bigPowChan = chan;
            end
        end
        fprintf('\n\nChannel with biggest oscillation power: %d', bigPowChan)

        baselineRawVec = baselineRawVec(bigPowChan, :);
        rawFilt = filtfilt(b, a, baselineRawVec);
        if peakPow == -1
            rawFilt = rawFilt .* (rawFilt < 0);
        end
        if peakPow == 1
            rawFilt = rawFilt .* (rawFilt > 0);
        end
        clear baselineRawVec
        rawPow = rawFilt.^2;
        clear rawFilt;
        if verbose; fprintf('\nConstructing envelope'); end
    %     [envUpper, ~] = envelope(rawPow, 10, 'peak');
        envUpper = getUpperEnvelope(rawPow, 'verbose', 1);
        clear rawPow;
        envMean = mean(envUpper);
        envStd = std(envUpper);
        clear envUpper;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('\n\nLooking at interest epoch(s)');
    if ~isnumeric(interestRawVec)
        if verbose; fprintf('\nReading in interest raw data'); end
        interestRawVec = readmda(char(interestRawVec));
    end
    interestRawVec = double(interestRawVec);
    
    if length(chanExamine) > 1
        fprintf('\nFinding channel with largest oscillation power');
        for chan = chanExamine
            fprintf('\nAnalysing channel: %d', chan)

            rawFilt = filtfilt(b, a, interestRawVec(chan, :));
            rawPow = rawFilt.^2;
            clear rawFilt;
            tempTotPow = sum(rawPow);
            if tempTotPow > bigPow
                bigPow = tempTotPow;
                bigPowChan = chan;
            end
        end
        fprintf('\n\nChannel with biggest oscillation power: %d', bigPowChan)
    else
        bigPowChan = chanExamine(1);
    end
    
    interestRawVec = interestRawVec(bigPowChan, :);
    rawFilt = filtfilt(b, a, interestRawVec);
    
    clear interestRawVec;
    if peakPow == -1
        rawFilt = rawFilt .* (rawFilt < 0);
    elseif peakPow == 1
        rawFilt = rawFilt .* (rawFilt > 0);
    elseif nargout > 3
        fullPeakPowDir = zeros(1, length(rawFilt));
        for idx = 2:(length(rawFilt) - 1)
            if (((rawFilt(idx) - rawFilt(idx - 1)) > 0) && ((rawFilt(idx + 1) - rawFilt(idx)) < 0))
                fullPeakPowDir(idx) = 1;
            end
            if (((rawFilt(idx) - rawFilt(idx - 1)) < 0) && ((rawFilt(idx + 1) - rawFilt(idx)) > 0))
                fullPeakPowDir(idx) = -1;
            end
        end
    end
    rawPow = rawFilt.^2;
%     clear rawFilt;
    if verbose; fprintf('\nConstructing envelope'); end
    envUpper = getUpperEnvelope(rawPow, 'verbose', 1);
    if isempty(envUpper); return; end
    
    if isempty(baselineRawVec)
        envMean = mean(envUpper);
        envStd = std(envUpper);
    end
    
    envBoolStrong = envUpper >= (envMean + (strongThresh * envStd));
    envBoolWeak = envUpper >= (envMean + (weakThresh * envStd));
    fprintf('\nLooking for oscillation events')
    while maxIdx < (length(envBoolStrong) - floor(fs / 1000)) %This ends the search a millisecond before the end of the recording so we don't run into edge cases
        maxIdx = maxIdx + 1;

        if envBoolStrong(maxIdx)
            minIdx = 0;
            minStartIdx = minIdx;
            minEndIdx = minIdx;
            peakIdx = maxIdx;
            while ~minStartIdx || ~minEndIdx
                minIdx = minIdx + 1;

                if envUpper(peakIdx) < envUpper(maxIdx + minIdx); peakIdx = maxIdx + minIdx; end %The index already starts at the first instance of being above the threshold, so the peak power can't be at an index before the first above threshold index

                if ~minStartIdx && (maxIdx == minIdx)
                    minStartIdx = minIdx - 1;
                end

                if ~minEndIdx && (maxIdx + minIdx > length(envBoolStrong))
                    minEndIdx = minIdx - 1;
                end

                if ~minStartIdx && ~envBoolWeak(maxIdx - minIdx); minStartIdx = minIdx; end
                if ~minEndIdx && ~envBoolWeak(maxIdx + minIdx); minEndIdx = minIdx; end
            end

            oscIdxVec = [oscIdxVec; [(maxIdx - minStartIdx) peakIdx (maxIdx + minEndIdx)]];

            maxIdx = (maxIdx + minEndIdx);

        end

    end
    
    if correctIdx
        fprintf('\n\nStarting oscillation curation step');
        oscIdxVec = idxCurate(envUpper, oscIdxVec, minLength, maxLength, minGap, 'verbose', verbose, 'develop', develop);        
    end
    
    if nargout > 1
        varargout{1} = envMean;
    end
    if nargout > 2
        varargout{2} = envStd;
    end
    if nargout > 3
        if peakPow ~= 0
            peakPowDir = ones(1, size(oscIdxVec, 1)) .* peakPow;
        else
            peakPowDir = zeros(1, size(oscIdxVec, 1));
            for idx = 1:size(oscIdxVec, 1)
                peakPowDir(idx) = fullPeakPowDir(oscIdxVec(idx, 2));
            end
        end
        varargout{3} = peakPowDir;
    end
    fprintf('\n');
end
