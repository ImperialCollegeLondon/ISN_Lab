clear all
close all
clc

%% Read Data
filename = 'eeg_1.mat';
load(filename);

%% Compute PSD
eeg_psd = eegSleepScorePSD(eeg);

%% Save structure
save(['sleep_psd','_',char(regexp(filename,'[0-9]','match')),'.mat'],'eeg_psd');

%% Peak Frequency
for ch = 1:eeg_psd.nc
    [eeg_psd.max_pxx{ch}, max_idx] = max(eeg_psd.psd{ch},[],2);
    eeg_psd.max_f{ch} = eeg_psd.freq(max_idx);
end
idx = linspace(1,length(eeg.scoring),length(eeg_psd.max_f{1}));
figure();
colorPlot(1:length(eeg_psd.max_f{1}),eeg_psd.max_f{1},eeg.scoring(floor(idx))+1,6);

%% Cross Correlation
for ch = 1:eeg_psd.nc
    for i = 1:max(eeg_psd.freq)
        for j = i:max(eeg_psd.freq)
            [xcf(i,j,:),lags,bounds(i,j,:)] = crosscorr(eeg_psd.psd{ch}(:,i),eeg_psd.psd{ch}(:,j));
        end
    end
end

[max_corr, idx] = max(xcf,[],3);
max_corr_lag = lags(idx);

%% Auto-Non-linear analysis (matlab)
psd_smooth = smooth(eeg_psd.psd{1}(1,:),3);
[pks, loc] = findpeaks(psd_smooth,eeg_psd.freq,'MinPeakProminence',0.00007,'MinPeakHeight',0.0002);
% plot(loc,pks,'o');

tol = 0.2;

for i = 1:length(loc)
    mul = floor(max(eeg_psd.freq)/loc(i));
    if mul >= 2
        % haromonic analysis
        for j = 2:mul
            harm = loc(closest(loc(i) * j,loc));
            if (harm - loc(i) * j) < tol
                harmonics.exist(i,j-1) = true;
                harmonics.loc(i,j-1) = harm;
            else
                harmonics.exist(i,j-1) = false;
                harmonics.loc(i,j-1) = 0;
            end
        end
    end
    
    % multiplexing
    for k = i+1:length(loc)
        sumPlex = loc(closest(loc(i) + loc(k),loc));
        diffPlex = loc(closest(abs(loc(i) - loc(k)),loc));
        
        if (sumPlex - (loc(i) + loc(k)) < tol)
            multiplexMatrix.exist(i,k) = true;
            multiplexMatrix.loc(i,k) = sumPlex;
        else
            multiplexMatrix.exist(i,k) = false;
            multiplexMatrix.loc(i,k) = 0;
        end
        if (abs(diffPlex - abs(loc(i) - loc(k))) < tol)
            multiplexMatrix.exist(k,i) = true;
            multiplexMatrix.loc(k,i) = diffPlex;
        else
            multiplexMatrix.exist(k,i) = false;
            multiplexMatrix.loc(k,i) = 0;
        end
        
    end
end

%% Cross non-linear analysis (matlab)
psd_smooth_1 = smooth(eeg_psd.psd{1}(19,:),3);
psd_smooth_2 = smooth(eeg_psd.psd{1}(20,:),3);
[p_1, loc_1] = findpeaks(psd_smooth_1,eeg_psd.freq,'MinPeakProminence',0.00007,'MinPeakHeight',0.0002);
[p_2, loc_2] = findpeaks(psd_smooth_2,eeg_psd.freq,'MinPeakProminence',0.00007,'MinPeakHeight',0.0002);
plot(loc_1,p_1,'o',loc_2,p_2,'o');

tol = 0.5;

for i = 1:length(loc_1)
    mul = floor(max(eeg_psd.freq)/loc_1(i));
    if mul >= 2
        % haromonic analysis
        for j = 2:mul
            harm = loc_2(closest(loc_1(i) * j,loc_2));
            if (harm - loc_1(i) * j) < tol
                harmonics.exist(i,j-1) = true;
                harmonics.loc(i,j-1) = harm;
            else
                harmonics.exist(i,j-1) = false;
                harmonics.loc(i,j-1) = 0;
            end
        end
    end
    
    % multiplexing
    for k = i+1:length(loc_1)
        sumPlex = loc_2(closest(loc_1(i) + loc_1(k),loc_2));
        diffPlex = loc_2(closest(abs(loc_1(i) - loc_1(k)),loc_2));
        
        if (sumPlex - (loc_1(i) + loc_1(k)) < tol)
            multiplexMatrix.exist(i,k) = true;
            multiplexMatrix.loc(i,k) = sumPlex;
        else
            multiplexMatrix.exist(i,k) = false;
            multiplexMatrix.loc(i,k) = 0;
        end
        if (diffPlex - abs(loc_1(i) - loc_1(k)) < tol)
            multiplexMatrix.exist(k,i) = true;
            multiplexMatrix.loc(k,i) = diffPlex;
        else
            multiplexMatrix.exist(k,i) = false;
            multiplexMatrix.loc(k,i) = 0;
        end
        
    end
end
%% Non-linear analysis (differential)
psd_smooth = smooth(eeg_psd.psd{1}(1,:),3);
df = eeg_psd.freq(2)-eeg_psd.freq(1);
dpsd = gradient(psd_smooth,df);
ddpsd = gradient(dpsd,df);
dddpsd = gradient(ddpsd,df);

loc = eeg_psd.freq(abs(dpsd) < 1E-1 & ddpsd < 0);
h = psd_smooth(abs(dpsd) < 1E-1 & ddpsd < 0 );

plot(eeg_psd.freq,psd_smooth);
hold on
plot(loc,h,'o');
hold off