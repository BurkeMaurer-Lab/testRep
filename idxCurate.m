%Version 1: Single use, non-modular

function [idxVec] = idxCurate(envUpper, idxVec, minLength, maxLength, minGap, varargin)

    %defaults
    varStrings = ["fs" "verbose" "develop"];
    fs = 24414;
    verbose = 0;
    develop = 0;
    for i = 1:2:length(varargin)
        if ~ismember(varargin{i}, varStrings) 
            fprintf('\n\nInput does not match allowable options.\nYou entered %s\nPlease try again.\n', string(varargin{i}))
            return;
        end
        eval([varargin{i} '=varargin{i + 1};']);
    end
    if develop; verbose = 1; end
    
    maxLengthIdx = floor(maxLength * (fs / 1000));
    minLengthIdx = floor(minLength * (fs / 1000));
    minGapIdx = floor(minGap * (fs / 1000));
    
    if develop
        preFig = figure;
        lengthVec = idxVec(:, 3) - idxVec(:, 1);
        lengthVec = lengthVec / (fs / 1000);
        histogram(lengthVec, 100)
        title('Pre Curate')
        hold on;
        vline(minLength, 'k');
        hold on;
        vline(maxLength, 'k');
        hold off;
    end
    
    idx = 1;
    numOsc = size(idxVec, 1);
    preCurateNum = numOsc;
    if verbose; fprintf('\nStarting to combine oscillation events'); end
    
    while idx < numOsc
        if develop && mod((numOsc - idx), 500) == 0; fprintf('\nNum Oscillations left to look at: %d', (numOsc - idx)); end
        if idxVec(idx + 1, 1) - idxVec(idx, 3) <= minGapIdx
            numOsc = numOsc - 1;
            newPowIdx = idxVec(idx, 2);
            if envUpper(idxVec(idx + 1, 2)) > envUpper(newPowIdx); newPowIdx = idxVec(idx + 1, 2); end
            idxVec(idx, :) = [idxVec(idx, 1) newPowIdx idxVec(idx + 1, 3)];
            idxVec(idx + 1, :) = [];
        else
            idx = idx + 1;
        end     
    end
    
    if develop
        postCombFig = figure;
        lengthVec = idxVec(:, 3) - idxVec(:, 1);
        lengthVec = lengthVec / (fs / 1000);
        histogram(lengthVec, 100)
        title('Post Combine')
        hold on;
        vline(minLength, 'k');
        hold on;
        vline(maxLength, 'k');
        hold off;
    end

    numOscPostComb = size(idxVec, 1);
    numOsc = numOscPostComb;
    idx = 1;

    if verbose; fprintf('\nStarting to cut out bad oscillations'); end
    while idx <= numOsc
        if develop && mod((numOsc - idx), 500) == 0; fprintf('\nNum Oscillations left to look at: %d', (numOsc - idx)); end
        if idxVec(idx, 3) - idxVec(idx, 1) < minLengthIdx || idxVec(idx, 3) - idxVec(idx, 1) > maxLengthIdx

            numOsc = numOsc - 1;
            if idx == 1
                idxVec = idxVec(2:end, :); 
            elseif idx == numOsc 
                idxVec = idxVec(1:(end - 1), :);
            else
                idxVec(idx, :) = [];
%                 = [idxVec(1:(idx - 1), :); idxVec((idx + 1):end, :)];
            end
        else
            idx = idx + 1;
        end
    end
    numOscPostCut = size(idxVec, 1);
    
    if develop
        postCutFig = figure;
        lengthVec = idxVec(:, 3) - idxVec(:, 1);
        lengthVec = lengthVec / (fs / 1000);
        histogram(lengthVec, 100)
        title('Post Combine and Cut')
        hold on;
        vline(minLength, 'k');
        hold on;
        vline(maxLength, 'k');
        hold off;
    end

    if verbose
        fprintf('\nStart Num Oscillations: %d', preCurateNum);
        fprintf('\nPost-Combine Num Oscillations: %d', numOscPostComb);
        fprintf('\nPost-Cut Num Oscillations: %d', numOscPostCut);
    end
        
    fprintf('\n')
end
