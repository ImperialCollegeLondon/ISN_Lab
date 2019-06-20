clear all
close all
clc

%% user input
sbjct_num = 3;
filename = '01-02-0003 PSG.edf';
filename_annotation = '01-02-0003 Base.edf';

% ch_to_use = [9, 16, 17, 22];    % Cz, O1, Pz, Fpz
ch_name = ["Cz","O1","Pz","Fpz"];
num_channels = length(ch_name);

% epoch size
ep_sz = 5; % if 0, follow annotation

%% Load eeg data
% load eeg file
[hdr, data] = edfread(filename);

%% Extract EEG data
ch_to_use = find(contains(hdr.label,ch_name)); % find eeg with matching name
if length(ch_to_use) < length(ch_name)
    warning('Not all channel exist')
    disp(hdr.label{ch_to_use});
end
eeg.data = data(ch_to_use,:)';  % eeg
eeg.nc = length(ch_to_use);
eeg.channels = cellfun(@(x) extractBetween(x,'EEG','CLE'),hdr.label(ch_to_use),'uni',false );
length_eeg = length(eeg.data);
Fs = hdr.frequency(ch_to_use);  % sampling rate
Fs = Fs(1);
dt = 1./Fs;
eeg.t = zeros(length_eeg, num_channels);
eeg.Fs = Fs;
for ch = 1:num_channels
    % channel label
    eeg.labels{ch} = hdr.label{ch_to_use(ch)} ;
end
% create t vector
eeg.t = [0:dt:dt*(length_eeg-1)];


%% extract anotation dta
annot = readEDFAnnotation(filename_annotation);

scoring = [zeros(length(annot.annotation),1)];

annot.annotation = regexp(annot.annotation,'[1-4]|(W)|(R)|(\?)','match');
annot.annotation = cellfun(@(x) regexprep(x,'W','0'),annot.annotation,'uni', false);
annot.annotation = cellfun(@(x) regexprep(x,'R','5'),annot.annotation,'uni', false);
annot.annotation = cellfun(@(x) regexprep(x,'\?','-1'),annot.annotation,'uni', false);

annot.annotation = vertcat(annot.annotation{:});
annot.annotation = cellfun(@str2num, annot.annotation,'uni', false);
annot.annotation = cell2mat(annot.annotation);

scoring = annot.annotation;

% put in stucture
eeg.scoring = scoring;

eeg.scoring_start = floor(annot.onset(1,1)*Fs)+1;
eeg.scoring_end = floor(annot.onset(end,1)*Fs+1 + annot.duration(end,1)*Fs);

% epoch of scoring
if ep_sz == 0
    eeg.ep_sz = floor(annot.duration(:,1) * Fs);
    eeg.ep_onset = floor(annot.onset(:,1)*Fs)+1;
else
    eeg.ep_onset = floor([eeg.scoring_start:ep_sz * Fs: eeg.scoring_end])';
    eeg.ep_sz = [floor(ep_sz * Fs)];
%     eeg.ep_onset = eeg.ep_onset(1:end-1);    
end

% %% Plot EEG with Score Colouring
% figure
% colorPlot(eeg.t,eeg.data(:,1),eeg.scoring+1,6);

%% save structure
save(['eeg','_',num2str(sbjct_num),'.mat'],'eeg');
fprintf('Done \n');