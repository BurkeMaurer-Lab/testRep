function [upEnvelope] = getUpperEnvelope(raw, varargin)

    upEnvelope = [];
    if size(raw, 1) ~= 1
        fprintf('\nMust be a one dimensional array. Please try again.')
        return;
    end
    
    %defaults
    varStrings = ["verbose" "fs"];
    verbose = 0;
    fs = 24414;
    for i = 1:2:length(varargin)
        if ~ismember(varargin{i}, varStrings) 
            fprintf('\n\nInput does not match allowable options.\nYou entered %s\nPlease try again.\n', string(varargin{i}))
            return;
        end
        eval([varargin{i} '=varargin{i + 1};']);
    end
    if verbose; fprintf('\n\nStarting envelope detection'); end
    
    upEnvelope = zeros(1, length(raw));
    prevBigV = raw(1);
    prevBigIdx = 1;
    for i = 1:length(raw)
        if verbose && (mod(i, (60 * fs)) == 0); fprintf('\nAt minute %d out of %d', (i / fs / 60), (length(raw) / fs / 60)); end
        if i == 1; upEnvelope(i) = raw(i); continue; end
        if (i == length(raw)) || (((raw(i) - raw(i - 1)) > 0) && ((raw(i + 1) - raw(i)) < 0))
            upEnvelope(i) = raw(i);
            curBigV = raw(i);
            fillVolt = linspace(prevBigV, curBigV, (i - prevBigIdx + 1));
            upEnvelope((prevBigIdx + 1):(i - 1)) = fillVolt(2:(end - 1));
            clear fillVolt;
            prevBigV = curBigV;
            prevBigIdx = i;                
        end
    end
    
    fprintf('\n')
end
