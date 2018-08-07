%Version 1: Only one file at a time

%Returns the indices of the start and stop of a sharp wave ripple along
%with the relative index of the maximum ripple power
%
%
%150 300 15 500 50
function [oscIdxVec] = getOscillationIndex(interestRawVec, baselineRawVec, highPass, lowPass, minLength, maxLength, minGap, varargin)
    
    oscIdxVec = [];
    bigPowEnvelope = [];
    bigPow = 0;
    bigPowChan = 0;

    %default variables
    varStrings = ["fs" "verbose" "develop" "swrMaxLength" "swrMinLength" "swrMinGap" "correctSWR" "passThresh" "minThresh" "chanExamin" "passType"];
    fs = 24414;
    correctIdx = 1; %If 1 the program will go through and combine any ripple events that have an intraripple time of less than or equal to the indicated gap length. With anything remaining, it'll then get rid of any SWR that are shorter than the minimum length or longer than the maximum length.
    verbose = 0;
    develop = 0;
    peakPow = -1; %-1 only looks at the negative deflecting waves for peak power, 0 looks at both, and 1 looks at the postive deflecting waves.
    passThresh = 5; %STD
    minThresh = 1; %IQR
    chanExamin = []; %Which channels should be looked at. An empty matrix means all of them.
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
    maxLengthIdx = maxLength * (fs / 1000);
    minLengthIdx = minLength * (fs / 1000);
    minGapIdx = minGap * (fs / 1000);
    
    if passType == 'bandpass'
        [b, a] = butter(2, [highPass lowPass] / (fs / 2), passType);
    elseif passType == 'lowpass'
        [b, a] = butter(2, [lowPass] / (fs / 2), passType);
    else
        [b, a] = butter(2, [highPass] / (fs / 2), passType);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FOR FUTURE DEVELOPMENT
%     if class(rawVec) == 'string'
%         numFiles = length(rawVec);
%     else
%         openedData = 1;
%         numFiles = size(rawVec, 3);
%     end

%     for fileIdx = 1:numFiles
%     if verbose; fprintf('\nAnalyzing file %d', fileIdx); end
%
%     clear raw;
%     if openedData
%         raw = rawVec(:, :, fileIdx);
%     else
%         raw = readmda(char(rawVec(fileIdx)));
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if class(baselineRawVec) ~= 'double'
        if verbose; fprintf('\nReading in behavior raw data'); end
        baselineRawVec = readmda(char(baselineRawVec));
    end

    if isempty(chanExamin); chanExamin = 1:size(baselineRawVec, 1); end

    fprintf('\n\nExamining baseline epoch to determine ripple thresholds.')
    fprintf('\nFinding channel with largest ripple power');
    for chan = chanExamin
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
    clear behaviorRawVec
    rawPow = rawFilt.^2;
    clear rawFilt;
    if verbose; fprintf('\nConstructing envelope'); end
%     [envUpper, ~] = envelope(rawPow, 10, 'peak');
    envUpper = getUpperEnvelope(rawPow, 'verbose', 1);
    clear rawPow;
    envMean = mean(envUpper);
    envSTD = std(envUpper);
    clear envUpper;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('\n\nLooking at rest epoch(s)');
    if class(interestRawVec) ~= 'double'
        if verbose; fprintf('\nReading in rest raw data'); end
        interestRawVec = readmda(char(interestRawVec));
    end
    fprintf('\nFinding channel with largest ripple power');
    for chan = chanExamin
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
    fprintf('\n\nChannel with biggest ripple power: %d', bigPowChan)
    
    interestRawVec = interestRawVec(bigPowChan, :);
    rawFilt = filtfilt(b, a, interestRawVec);
    
    clear restRawVec;
    if peakPow == -1
        rawFilt = rawFilt .* (rawFilt < 0);
    end
    if peakPow == 1
        rawFilt = rawFilt .* (rawFilt > 0);
    end
    rawPow = rawFilt.^2;
%     clear rawFilt;
    if verbose; fprintf('\nConstructing envelope'); end
    envUpper = getUpperEnvelope(rawPow, 'verbose', 1);
    if isempty(envUpper); return; end
    
    envBoolMax = envUpper >= (envMean + (5 * envSTD));
    envBoolMin = envUpper >= (envMean + envSTD);
    fprintf('\nLooking for SWR events')
    while maxIdx < (length(envBoolMax) - floor(fs / 1000)) %This ends the search a millisecond before the end of the recording so we don't run into edge cases
        maxIdx = maxIdx + 1;

        if envBoolMax(maxIdx)
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

                if ~minEndIdx && (maxIdx + minIdx > length(envBoolMax))
                    minEndIdx = minIdx - 1;
                end

                if ~minStartIdx && ~envBoolMin(maxIdx - minIdx); minStartIdx = minIdx; end
                if ~minEndIdx && ~envBoolMin(maxIdx + minIdx); minEndIdx = minIdx; end
            end

            oscIdxVec = [oscIdxVec; [(maxIdx - minStartIdx) peakIdx (maxIdx + minEndIdx)]];

            maxIdx = (maxIdx + minEndIdx);

        end

    end
    
    if correctIdx
        fprintf('\n\nStarting ripple curation step');
        oscIdxVec = idxCurate(oscIdxVec, minLength, maxLength, minGap, 'verbose', verbose, 'develop', develop);        
    end
        
    fprintf('\n');
end
