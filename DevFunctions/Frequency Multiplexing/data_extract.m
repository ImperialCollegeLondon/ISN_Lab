function data = data_extract(filename, annotation)

% Reading Header
hdr = edfread(filename);
% Predefined index for EEG channels
% target = [2,3,5,8,9,10,11,13,14,15,16,17,18,19,20,22,23,24,25,26];
target = find(contains(hdr.label,'EEG'));

% Reading in the EDF
[hdr, data] = edfread(filename,'targetSignals',target);

if size(data,1) < 20
    fprintf('Not all data channels are present!');
    disp(hdr.label);
end

% Export channels into text file
for i = 1:size(data,1)
    fprintf('Writing channel %s',hdr.label{target(i)});
    fileID = fopen(strcat(extractBefore(filename,'.edf'),'_',hdr.label{target(i)},'.txt'),'w');
    fprintf(fileID,'Fs: %.4f\n',hdr.frequency);
    fprintf(fileID,'%6s %12s\n','time (ms)',strcat('amplitude (',hdr.units{target(i)},')'));
    fprintf(fileID,'%6d %12.8f\n',[1:length(data(i,:))],data(i,:));
    fclose(fileID);
end

% Not Ready (Requires pre-extracted annotations)
% if nargin == 2
%     % Annotation-based extraction
%     scoring_f = annotation;
%     
%     Fs = hdr.frequency;
%     
%     % import recorded start and duration of spindles
%     delimiter = ',';
%     startRow = 2;
%     formatSpec = '%f%f%*s%[^\n\r]';
%     fileID = fopen(scoring_f,'r');
%     scoring = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     fclose(fileID);
%     
%     classification = zeros(size(data)); % Binary vector for spindle/baseline seperation
%     
%     % Set start and finish point of spindles
%     start = floor(scoring{1,1}*Fs);
%     duration = floor(scoring{1,2}*Fs);
%     finish = start+duration;
%     
%     spindle = cell(1,length(start));
%     baseline = cell(1,length(start));
%     
%     % Classify the time series and pickout spindle segments. Baseline segment
%     % is randomly assigned as segments before the spindle.
%     for i = 1:length(start)
%         classification(start(i):finish(i)) = 1;
%         spindle{i} = data(:,start(i):finish(i));
%         
%         a = finish(max(1,i-1))+1;
%         
%         r = a + round((start(i)-1-a).*rand);
%         while classification(r:r+floor(Fs/2)-1) == 1
%             r = a + round((start(i)-1-a).*rand);
%         end
%         baseline{i} = data(:,r:r+floor(Fs/2)-1);
%     end
%     
%     % Truncate spindle data to 50 samples (centered)
%     baseline(cellfun(@(x) size(x,2),spindle) < floor(Fs/2)) = [];
%     spindle(cellfun(@(x) size(x,2),spindle) < floor(Fs/2)) = [];
%     spindle = cellfun(@(x) x(:,floor(end/2)-floor(Fs/4-1):floor(end/2)+floor(Fs/4)),spindle,'UniformOutput',false);
%     
%     clearvars -except baseline spindle Fs
%     
% end

end