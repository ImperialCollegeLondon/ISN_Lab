function [spindle, baseline] = data_extract_ready(data, scoring_f, Fs)

% import recorded start and duration of spindles
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%*s%[^\n\r]';
fileID = fopen(scoring_f,'r');
scoring = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

classification = zeros(size(data)); % Binary vector for spindle/baseline seperation

% Set start and finish point of spindles
start = floor(scoring{1,1}*Fs);
duration = floor(scoring{1,2}*Fs);
finish = start+duration;

spindle = cell(1,length(start));
baseline = cell(1,length(start));

% Classify the time series and pickout spindle segments. Baseline segment
% is randomly assigned as segments before the spindle.
for i = 1:length(start)
    a = -1;
    classification(start(i):finish(i)) = 1;
    spindle{i} = data(start(i):finish(i));
    
    if i == 1
        a = 1;
    else
        a = finish(i-1)+1;
    end
    
    r = a + round((start(i)-1-a).*rand);
    while classification(r:r+floor(Fs/2)-1) == 1
        r = a + round((start(i)-1-a).*rand);
    end
    baseline{i} = data(r:r+floor(Fs/2)-1);
end

% Truncate spindle data to 50 samples (centered)

clearvars -except baseline spindle Fs