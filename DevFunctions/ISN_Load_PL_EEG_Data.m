function [dataArray] = ISN_Load_PL_EEG_Data(filename)
%ISN_Load_PL_EEG_Data Loads the data from the output .csv file from our PL_EEG system
%   [dataArray] = ISN_Load_PL_EEG_Data(filename); Current Firmware: 190307
%   This function takes into account the three different firmware used to
%   create EEG files so far. The output variable always has the following
%   order: (1) Abs time, (2) Mode, (3) DT, (4) EEG, (5) Stim_out

% Makes sure the filename provided is for a .csv file
if numel(strfind(filename,'.csv')) ==1
    tmp_filename = filename;
else
    tmp_filename = [filename '.csv'];
end

% Determines the Firmware used to create the file then imports the data to
% dataArray
delimiter = ',';
fileID = fopen(tmp_filename,'r');
topline = fgets(fileID);
numCols = numel(strfind(topline,delimiter)) + 1;
if numCols == 7 % FIRMWARE 181210
    disp ('Firmware version: 181210')
    formatSpec = '%f%f%f%f%f%*s%*s%[^\n\r]';
    dataArray2 = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
elseif numCols == 11 % FIRMWARE 190205 (beta)
    disp ('Firmware version: 190205')
    formatSpec = '%f%*s%f%f%f%f%*s%*s%*s%*s%*s%[^\n\r]';
    dataArray2 = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
elseif numCols == 9 %FIRMWARE 190307
    disp ('Firmware version: 190307')
    formatSpec = '%*s%f%*s%f%f%f%f%*s%*s%*s%*s%*s%[^\n\r]';
    dataArray2 = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
end
fclose(fileID);
dataArray = [dataArray2{1:end-1}];

% Firmware version 181210 is in a different order to the others so needs to
% be re-ordered
if numCols == 7
    tmp = [dataArray(:,1), dataArray(:,5), dataArray(:,4), dataArray(:,2), dataArray(:,3)];
    dataArray = tmp;
    % Firmware version 190307 now contains a log - @echo, these lines do
    % not contain data and have been removed
elseif numCols == 9
    for n = size(dataArray,2):-1:1
        tmp(:,n) = isnan(dataArray(:,n));
    end
    tmp2=sum(tmp,2);
    tmp3 = tmp2==0;
    dataArray = dataArray(tmp3,:);
end

end

