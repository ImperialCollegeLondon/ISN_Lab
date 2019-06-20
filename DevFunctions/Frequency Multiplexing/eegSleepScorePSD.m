function psd = eegSleepScorePSD(eeg)
%eegPSD plot PSD and returns eeg struct with PSD
%
%% USER INPUT

% if ~isempty(varargin)
%     for ii = 1:2:numel(varargin)
%         switch lower(varargin{ii})
%             case 'psd_subject_level_analysis'
%                 psd_subject_level_analysis = varargin{ii+1};
%             case 'targetsignals'
%                 targetSignals = varargin{ii+1};
%             otherwise
%                 error('eegPSD: Unrecognized parameter-value pair specified.')
%         end
%     end
% end

%%% experiment design (conditions, locations, periods)
dscrp_loc = eeg.channels;
num_locations = eeg.nc;

% normalized to total power
normalized_psd_to_total_power = 1;

%% COMPUTE EEG POWER SPECTRAL DENSITY (PSD)

% extract subject data
y_eeg = eeg.data;
Fs = eeg.Fs;
num_epochs = length(eeg.ep_onset);
if length(eeg.ep_sz) > 2
    epoch_idx = arrayfun(@(x,y) x:x+y, eeg.ep_onset,eeg.ep_sz,'uni',false);
else
    for i = 1:length(eeg.ep_onset)
        epoch_idx{i} = [eeg.ep_onset(i):eeg.ep_onset(i) + eeg.ep_sz(1)];
    end
end
% excld_epoch = (eeg.scoring == -1);

%limit freq span
max_pxx_f = 30;

% Temporary variables
tmp_pxx = [];
pxx = [];
f = [];
tmp_f = [];

nfft = 2^(nextpow2(eeg.ep_sz(1)));

for ch=1:eeg.nc
    % compute power density spectrum at epochs
    fprintf('Evaluating channel %d\n',ch);
    tic
    for i=1:num_epochs
        % extract epoch data_all_sbjcts
        yy = y_eeg(epoch_idx{i},ch);
        % compute multitaper power spectral density (dB/Hz) in epoch
        [tmp_pxx,tmp_f] = pmtm(yy,3,nfft, Fs); % length of yy may vary per epoch
        % set idx for the 30 Hz mark
        idx = round(max_pxx_f/(Fs/nfft),-1);
        pxx{ch}(i,:) = tmp_pxx(1:idx);
        % normalized to total power
        if(normalized_psd_to_total_power)
            pxx{ch}(i,:) = pxx{ch}(i,:)/sum(pxx{ch}(i,:));
        end
        % set epoch with artifact to nan
%         if(excld_epoch(epoch_idx{i}))
%             pxx{ch}(i,:) = NaN;
%         end
    end
    toc
end

f = tmp_f(1:idx);

%% Store into variable
psd.nc = eeg.nc;
psd.nepc = num_epochs;
psd.channels = eeg.channels;
psd.freq = f;
psd.psd = pxx;
psd.scoring = eeg.scoring;

end

