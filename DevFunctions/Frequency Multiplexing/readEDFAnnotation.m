function annotation = readEDFAnnotation(filename)
%readEDFAnnotation Reads Annotation from EDF file and store them as vector
%   filename : Input file name
%
%   annotation : Outputs a struct that include onset, duration and
%   its respective annotation.
%   For spindle annotation, 0 is baseline, 1 is spindle and -1 is no
%   annotation
%   For sleep score: 0 = wake
%                    1 = N1 NREM/light sleep
%                    2 = Sigma/N2 NREM/Spindle Sleep
%                    3 = N3 NREM/Deep Sleep/Slow-wave Sleep
%                    4 = N4 NREM/Slow-wave Sleep
%                    5 = REM
%                    -1 = No Annotations
%
%   Author: Lok Chun Fan
%   Date: 22/2/2019

% Open EDF file
[fid,msg] = fopen(filename,'r','l','UTF-8');
if fid == -1
    error(msg)
end

fseek(fid,236,-1);
records = str2double(fread(fid,8,'*char')');
fseek(fid,472,-1);
ns = str2double(fread(fid,8,'*char')');

annotation.onset = zeros(records,1);
annotation.duration = zeros(records,1);
annotation.annotation = cell(records,1);

% Reading annotation samples
fseek(fid,512,-1);
onset = 0;
duration = 0;
onsetTime = 0;
durationTime = 0;
tmp = zeros(1,20);
i = 1;

for j = 1: records
    emptyByte = 0;
    
    for k = 1:ns*2
        curChar = fread(fid, 1, '*char');
        if curChar == '+'
            onset = 1;
            i = 1;
            continue
        end
        
        if onset == 1
            if curChar == 20
                onset = 0;
                tmp = zeros(1,20);
                i = 1;
            elseif curChar == 21
                onset = 0;
                duration = 1;
                onsetTime = str2double(char(tmp));
                tmp = zeros(1,20);
                i = 1;
            else
                if i > 20
                    error('Exceeded buffer size')
                end
                tmp(i) = curChar;
                i = i + 1;
            end
        elseif duration == 1
            if curChar == 20
                duration = 0;
                durationTime = str2double(char(tmp));
                tmp = zeros(1,20);
                i = 1;
            else
                if i > 20
                    error('Exceeded buffer size')
                end
                tmp(i) = curChar;
                i = i + 1;
            end
        elseif curChar == 0
            emptyByte = emptyByte + 1;
            if emptyByte == 2
                break
            end
        else
            if curChar == 20 
                continue 
            end
            annotation.annotation{j} = strcat(annotation.annotation{j,1},curChar);
        end
    end
    
    if ~isempty(annotation.annotation{j,1})
        annotation.onset(j) = onsetTime;
        annotation.duration(j) = durationTime;
        onsetTime = 0;
        durationTime = 0;
    end
    
    % Jump to next block of record
    fseek(fid,512+(j)*ns*2,-1);
end

fclose(fid);

annotation.annotation = annotation.annotation(~cellfun('isempty',annotation.annotation));
annotation.onset = annotation.onset(annotation.onset ~= 0);
annotation.duration = annotation.duration(annotation.duration ~= 0);

end