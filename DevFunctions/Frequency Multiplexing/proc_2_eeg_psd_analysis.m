clear all
close all
clc

%% Read Data
filename = 'eeg_1.mat';
load(filename);

%% Modify Epoch Size
ep_sz = 5;

eeg.ep_onset = floor([eeg.scoring_start:ep_sz * eeg.Fs: eeg.scoring_end])';
eeg.ep_sz = [floor(ep_sz * eeg.Fs)];

%% Compute PSD
cmp_psd = 1;
eeg_psd = struct();
if(cmp_psd)
    eeg_psd = eegSleepScorePSD(eeg);
    
    % Save structure
    save(['sleep_psd','_',char(regexp(filename,'[0-9]','match')),'.mat'],'eeg_psd');
else
    % load PSD strucure
    load(['sleep_psd','_',char(regexp(filename,'[0-9]','match')),'.mat']);
end

%% Channel and Epoch Settings
ch_start = 1;
ch_end = eeg_psd.nc;

epch_start = 1;
epch_end = eeg_psd.nepc;

if ch_end > eeg_psd.nc
    disp("channel range larger than number of channels, using max channel instead")
    ch_end = eeg_psd.nc;
elseif ch_end < ch_start
    error("channel range cannot be negative")
end

if epch_end > eeg_psd.nepc
    disp("epoch range larger than number of epochs, using max epoch instead")
    epch_end = eeg_psd.nepc;
elseif epch_end < epch_start
    error("epoch range cannot be negative")
end

%% Auto-Non-linear analysis (matlab)
close all
tol = 0.03; %TODO: percentage tolerance

max_harmonic = 3;

ch_range = ch_end-ch_start + 1;
epch_range = epch_end - epch_start + 1;

eeg_psd.pks_amp = cell(ch_range,epch_range);
eeg_psd.pks_freq = cell(ch_range,epch_range);
eeg_psd.npks = zeros(ch_range,epch_range);

eeg_multiplex.pks_amp = cell(ch_range,epch_range);
eeg_multiplex.pks_freq = cell(ch_range,epch_range);
eeg_multiplex.triplet_count = cell(ch_range,epch_range);
eeg_multiplex.triplet_component = cell(ch_range,epch_range);
% eeg_multiplex.diff_count = cell(ch_range,epch_range);
% eeg_multiplex.diff_contribution = cell(ch_range,epch_range);
eeg_multiplex.harmonic_count = cell(ch_range,epch_range);
eeg_multiplex.harmonic_component = cell(ch_range,epch_range);

tic
% main loop (put instructions to loop through channels and epochs here)
for ch = ch_start:ch_end
    fprintf("Evaluating channel %d\n",ch)
    for epch = epch_start:epch_end
        
        % Signal smoothing and Peak finding
        num_psd_smooth_points = 10;
        psd_smooth = smooth(eeg_psd.psd{ch}(epch,:),num_psd_smooth_points);
        %         [eeg_psd.pks_amp{ch,epch}, eeg_psd.pks_freq{ch,epch}] = findpeaks(psd_smooth, eeg_psd.freq, 'MinPeakProminence',0.00007,'MinPeakHeight',0.0002);
        [eeg_psd.pks_amp{ch,epch}, eeg_psd.pks_freq{ch,epch}] = findpeaks(psd_smooth, eeg_psd.freq,'threshold',0.0000001, 'MinPeakProminence',0.0001,'MinPeakHeight',0.0001);
        
        % compute num of peaks
        eeg_psd.npks(ch,epch) = length(eeg_psd.pks_freq{ch,epch});
        
        [triplet_count, triplet_component, harmonic_count, harmonic_component] = multiplex_find(eeg_psd.pks_freq{ch,epch},eeg_psd.pks_freq{ch,epch}, tol, max(eeg_psd.freq), 3);
        
        
        % saving data into structure
        eeg_multiplex.triplet_count{ch,epch} = [eeg_psd.pks_freq{ch,epch} triplet_count'];
        % eeg_multiplex.diff_count{ch,epch} = [eeg_psd.pks_freq{ch,epch} pks_diff_count'];
        eeg_multiplex.triplet_component{ch,epch} = [num2cell(eeg_psd.pks_freq{ch,epch}) triplet_component];
        % eeg_multiplex.diff_contribution{ch,epch} = [num2cell(eeg_psd.pks_freq{ch,epch}) pks_diff_contribution];
        eeg_multiplex.harmonic_count{ch,epch} = [eeg_psd.pks_freq{ch,epch} harmonic_count'];
        eeg_multiplex.harmonic_component{ch,epch} = [num2cell(eeg_psd.pks_freq{ch,epch}) harmonic_component];
        
    end
end

toc

% saving peaks information in multiplex data structure
eeg_multiplex.pks_amp = eeg_psd.pks_amp;
eeg_multiplex.pks_freq = eeg_psd.pks_freq;

% saving sleep score data into multiplex data structure
eeg_multiplex.sleep_score = eeg_psd.scoring;

% saving other useful information
eeg_multiplex.nc = eeg_psd.nc;
eeg_multiplex.nepc = eeg_psd.nepc;
eeg_multiplex.npks = eeg_psd.npks;
eeg_multiplex.channels = eeg_psd.channels;

%% Duo-epoch Multiplexing Analysis

eeg_multiplex = duoEpochMultiplex(eeg_multiplex,eeg_multiplex.nc,eeg_multiplex.nepc,eeg_psd.pks_freq, max(eeg_psd.freq));


%% Cleanup variables
clearvars -except eeg* ch_* epch_* filename

%% Save multiplex data structure
save(['sleep_multiplex','_',char(regexp(filename,'[0-9]','match')),'.mat'],'eeg_multiplex');
save(['sleep_psd','_',char(regexp(filename,'[0-9]','match')),'.mat'],'eeg_psd');

%% Peaks visualisation
ch = 1;
epch = 1;
num_psd_smooth_points = 10;

psd_smooth = smooth(eeg_psd.psd{ch}(epch,:),num_psd_smooth_points);

figure
hold on
plot(eeg_psd.freq,psd_smooth,'k');
plot(eeg_multiplex.pks_freq{ch,epch}, eeg_multiplex.pks_amp{ch,epch},'ro','MarkerSize',5)
for i=1:length(eeg_multiplex.pks_freq{ch,epch})
    text(eeg_multiplex.pks_freq{ch,epch}(i), eeg_multiplex.pks_amp{ch,epch}(i)+0.0005, ['(',num2str(i),', ',num2str(eeg_multiplex.pks_freq{ch,epch}(i)),')'],'FontSize', 12)
end
hold off
set(gca, 'YScale', 'log')