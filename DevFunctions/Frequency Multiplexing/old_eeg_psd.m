clc
close all
clear all

%% NOTES
% ADC input range: +/- 2048
% DAC output range: 1-4096 but after loading it is centered to +/- 2048

%% USER INPUT
%%% subjects
sbjcts = [5:10,12:19];
% sbjcts = [5:6];

%%% ECHT setting
phase_lock_cntr_freq  = 5;

%%% experiment design (conditions, locations, periods)
dscrp_cond = {'peak', 'trough', 'sham'};
dscrp_loc = {'Oz', 'Fz'};
dscrp_prd = {'Pre', 'EO', 'EC'};
num_conditions = 3;
num_locations = 2;
num_periods = 3;

%%% PSD subject-level
psd_subject_level_analysis = 0;
% normalized to total power
normalized_psd_to_total_power = 0;
% z-transform
z_transform_psd = 1;    % 0,no; 1, yes
% sham-substra
substract_sham_psd = 1; % 0,no; 1, yes
% mean vs median
use_sbjct_median = 1;
% summary structure
save_all_sbjcts_psd_structure = 0;

%%% PSD group-level
psd_group_level_analysis = 1;
load_subject_level_structure = 1;
save_grp_psd_structure = 0;

%%% folder path
currentFolder = pwd;

%%% plotting setting
plot_sbjct_psd = 0;
plot_grp_psd = 1;
saving_plot = 1;


%% PLOT PROPERTIES
PlotLineWidth = 2.5;
PlotFontSize = 22;
PlotAxesLineWidth = 2;
PlotMarkerSize = 10;
set(0,'defaultLineLineWidth',PlotLineWidth);   % set the default line width to lw
set(0,'defaultLineMarkerSize',PlotMarkerSize); % set the default line marker size to msz
PlotPrintResolution = '-r300';
ColorLine = {[0.65, 0.1, 0.2],[0, 0.35, 0.65]};


%% %%%%%%%%%%% PSD SUBJECT-LEVEL %%%%%%%%%%%%%%
if(psd_subject_level_analysis)
    disp('%%%%%%%%%%%%% PSD SUBJECT-LEVEL %%%%%%%%%%%%%%');
    
    %% LOAD SUBJECT STRUCTURE
    DataPath = [currentFolder,'\1 - subject_level'];
    AnalysisPath = [currentFolder,'\1 - subject_level'];
    
    cd(DataPath);
    load('data_all_sbjcts.mat');
    
    cd(AnalysisPath);
    
    num_sbjcts = length(sbjcts);
    for n=1:num_sbjcts
        sbjct_num = sbjcts(n);
        disp(['* ',num2str(n),'. subject number = ', num2str(sbjct_num)])
        
        %% COMPUTE EEG POWER SPECTRAL DENSITY (PSD)
        tmp_pxx = [];
        pxx_baseline_mean = [];
        pxx = [];
        f = [];
        tmp_f = [];
        
        % extract subject data
        y_eeg = data_all_sbjcts.eeg.sbjct{n};
        Fs = data_all_sbjcts.eeg_Fs.sbjct{n};
        epochs_indxs = data_all_sbjcts.epochs_indxs.sbjct{n};
        num_epochs = data_all_sbjcts.num_epochs.sbjct{n}; 
        incl_epochs = data_all_sbjcts.incl_epochs.sbjct{n};  
        period_epochs = data_all_sbjcts.period_epochs.sbjct{n};
        
        
        for loc=1:num_locations
            for cond=1:num_conditions
                % compute power density spectrum at epochs
                for i=1:num_epochs.loc{loc}.cond{cond}
                    % extract epoch data_all_sbjcts
                    yy = y_eeg.loc{loc}.cond{cond}(epochs_indxs.loc{loc}.cond{cond}(i,1) : epochs_indxs.loc{loc}.cond{cond}(i,2));
                    % compute multitaper power spectral density (dB/Hz) in epoch
                    [tmp_pxx.loc{loc}.cond{cond}(i,:),tmp_f.loc{loc}.cond{cond}(i,:)] = pmtm(yy,3,length(yy), Fs.loc{loc}.cond{cond});
                    % limit freq span
                    max_pxx_f = 30;
                    pxx.loc{loc}.cond{cond}(i,:) = tmp_pxx.loc{loc}.cond{cond}(i,1:max_pxx_f);
                    f.loc{loc}.cond{cond}(i,:) = round(tmp_f.loc{loc}.cond{cond}(i,1:max_pxx_f), 1);
                    % normalized to total power
                    if(normalized_psd_to_total_power)
                        pxx.loc{loc}.cond{cond}(i,:) = pxx.loc{loc}.cond{cond}(i,:)/sum(pxx.loc{loc}.cond{cond}(i,:));
                    end
                    % set epoch with artifact to nan
                    if(~incl_epochs.loc{loc}.cond{cond}(i))
                        pxx.loc{loc}.cond{cond}(i,:) = NaN;
                    end
                end
            end
        end
        
        %% Z-TRANSFORM EEG PSD TO BASELINE
        if(z_transform_psd)
            for loc=1:num_locations
                for cond=1:num_conditions
                    pxx_baseline_mean.loc{loc}.cond{cond} = nanmean(pxx.loc{loc}.cond{cond}(period_epochs.loc{loc}.cond{cond}{1},:));
                    pxx_baseline_std.loc{loc}.cond{cond} = nanstd(pxx.loc{loc}.cond{cond}(period_epochs.loc{loc}.cond{cond}{1},:));
                    for i=1:num_epochs.loc{loc}.cond{cond}
                        pxx.loc{loc}.cond{cond}(i,:) = (pxx.loc{loc}.cond{cond}(i,:)-pxx_baseline_mean.loc{loc}.cond{cond})./pxx_baseline_std.loc{loc}.cond{cond};
                    end
                end
            end
        end
        
        %% COMPUTE PSD STATS IN PERIODS
        pxx_period = [];
        pxx_f_period = [];
        pxx_period_std = [];
        pxx_f_period_std = [];
        period_num_epochs_incl = [];
        for loc=1:num_locations
            for cond=1:num_conditions
                for prd=1:num_periods
                    pxx_period.loc{loc}.period{prd}.cond{cond} = nanmean(pxx.loc{loc}.cond{cond}(period_epochs.loc{loc}.cond{cond}{prd},:));
                    if(use_sbjct_median)
                        pxx_period.loc{loc}.period{prd}.cond{cond} = nanmedian(pxx.loc{loc}.cond{cond}(period_epochs.loc{loc}.cond{cond}{prd},:));
                    end
                    pxx_period_std.loc{loc}.period{prd}.cond{cond} = nanstd(pxx.loc{loc}.cond{cond}(period_epochs.loc{loc}.cond{cond}{prd},:));
                    period_num_epochs_incl.loc{loc}.period{prd}.cond{cond} = sum(~isnan(pxx.loc{loc}.cond{cond}(period_epochs.loc{loc}.cond{cond}{prd},1)));
                    pxx_f_period.loc{loc}.period{prd}.cond{cond} = nanmean(f.loc{loc}.cond{cond}(period_epochs.loc{loc}.cond{cond}{prd},:));
                    pxx_f_period_std.loc{loc}.period{prd}.cond{cond} = nanstd(f.loc{loc}.cond{cond}(period_epochs.loc{loc}.cond{cond}{prd},:));
                end
            end
        end
        
        
        %% SUBSTRACT SHAM PSD (per period)
        if(substract_sham_psd)
            for loc=1:num_locations
                for prd=1:num_periods
                    sham_mean = pxx_period.loc{loc}.period{prd}.cond{3};
                    for cond=1:num_conditions
                        pxx_period.loc{loc}.period{prd}.cond{cond} = ...
                            pxx_period.loc{loc}.period{prd}.cond{cond} - sham_mean;
                    end
                end
            end
        end
        
        %% PLOT SBJCT PSD
        if(plot_sbjct_psd)
            for loc=1:num_locations
                for prd=2:3 % EO and EC (basline is not plotted)
                    h1 = figure;
                    title = ['pwr_spctrm_sbj_',num2str(sbjct_num),'_', dscrp_loc{loc},'_', dscrp_prd{prd}];
                    set(h1,'name',title,'numbertitle','off');
                    set(gca, 'FontSize', PlotFontSize, 'LineWidth',PlotAxesLineWidth); %<- Set properties
                    for cond=1:2
                        hold on
                        % extract vectors
                        v_f = pxx_f_period.loc{loc}.period{prd}.cond{cond};
                        v_mean = pxx_period.loc{loc}.period{prd}.cond{cond};
                        %                     v_err = pxx_period_std.loc{loc}.period{prd}.cond{cond}; % std
                        v_err = pxx_period_std.loc{loc}.period{prd}.cond{cond}/sqrt(period_num_epochs_incl.loc{loc}.period{prd}.cond{cond}); % sem
                        % plot mean
                        plot(v_f, v_mean ,'color',ColorLine{cond},'LineWidth',2,'MarkerFaceColor',ColorLine{cond});
                        % plot std
                        errorbar(v_f, v_mean, v_err,'color',ColorLine{cond},'LineWidth',2);
                        hold off
                        xlim([0 25]);
                        %                     ylim([0 pxx_max+0.001]);
                        set(gca,'TickDir','out');
                        set(gca,'ticklength',3*get(gca,'ticklength'))
                        xlabel({'Frequency (Hz)'})
                        ylabel({'PSD (z-score)'})
                        if (saving_plot)
                            set(gcf,'PaperUnits','inches');
                            set(gcf,'PaperPosition',[0 0 11 7]);
                            set(gcf,'PaperPositionMode','manual')
                            name = [title,'.fig'];
                            saveas(h1,name);
                            name = [title];
                            print(name,'-dtiff',PlotPrintResolution);
                        end
                    end
                end
            end
        end
        
        
        %% CREATE ALL-SUBJECTS SUMMARY STRUCTURE
        % put in structure
        for loc=1:num_locations
            for prd=1:num_periods
                for cond=1:num_conditions
                    % PSD in epochs
                    data_all_sbjcts_psd.pxx.loc{loc}.cond{cond}.sbjct{n} = pxx.loc{loc}.cond{cond};
                    % PSD in periods
                    data_all_sbjcts_psd.pxx_period.loc{loc}.period{prd}.cond{cond}(n,:) = pxx_period.loc{loc}.period{prd}.cond{cond};
                    data_all_sbjcts_psd.pxx_period_std.loc{loc}.period{prd}.cond{cond}(n,:) = pxx_period_std.loc{loc}.period{prd}.cond{cond};
                    % PSD freq vector in periods
                    data_all_sbjcts_psd.pxx_f_period.loc{loc}.period{prd}.cond{cond}(n,:) = pxx_f_period.loc{loc}.period{prd}.cond{cond};
                    data_all_sbjcts_psd.pxx_f_period_std.loc{loc}.period{prd}.cond{cond}(n,:) = pxx_f_period_std.loc{loc}.period{prd}.cond{cond};
                end
            end
        end
        
        
    end % end 'n=1:num_sbjcts'
    
    %% SAVE STRUCTURE
    if(save_all_sbjcts_psd_structure)
        save(['data_all_sbjcts_psd','.mat'],'data_all_sbjcts_psd');
    end
    
end % end 'subect_level_analysis'


%% %%%%%%%%%%% PSD GROUP-LEVEL %%%%%%%%%%%%%%
if(psd_group_level_analysis)
    disp('%%%%%%%%%%%%% PSD GROUP-LEVEL %%%%%%%%%%%%%%');
    %% LOAD PSD SUBJECT STRUCTURE
    DataPath = [currentFolder,'\1 - subject_level'];
    AnalysisPath = [currentFolder,'\2 - group_level'];
    
    if(load_subject_level_structure)
        cd(DataPath);
        load('data_all_sbjcts_psd.mat');
    end
    
    cd(AnalysisPath);
    
    %% STATS
    sum_array_grp_psd = [];
    data_grp_psd = [];
    i = 1;
    for loc=1:num_locations
        for prd=1:num_periods
            for cond=1:num_conditions
                % psd
                tmp = data_all_sbjcts_psd.pxx_period.loc{loc}.period{prd}.cond{cond};
                data_grp_psd.stats.mean.pxx_period.loc{loc}.period{prd}.cond{cond} = mean(tmp);
                data_grp_psd.stats.median.pxx_period.loc{loc}.period{prd}.cond{cond} = median(tmp);
                data_grp_psd.stats.std.pxx_period.loc{loc}.period{prd}.cond{cond} = std(tmp);
                data_grp_psd.stats.n.pxx_period.loc{loc}.period{prd}.cond{cond} = size(tmp,1);
                % f
                tmp = data_all_sbjcts_psd.pxx_f_period.loc{loc}.period{prd}.cond{cond};
                data_grp_psd.stats.mean.pxx_f_period.loc{loc}.period{prd}.cond{cond} = mean(tmp);
                data_grp_psd.stats.std.pxx_f_period.loc{loc}.period{prd}.cond{cond} = std(tmp);
                
                % put in a summary array
                sum_array_grp_psd.pxx_period.row_dscrp{i,1} = [dscrp_loc{loc},'_', dscrp_prd{prd},'_',dscrp_cond{cond}];
                sum_array_grp_psd.pxx_period.mean(i,:) = data_grp_psd.stats.mean.pxx_period.loc{loc}.period{prd}.cond{cond};
                sum_array_grp_psd.pxx_period.median(i,:) = data_grp_psd.stats.median.pxx_period.loc{loc}.period{prd}.cond{cond};
                sum_array_grp_psd.pxx_period.std(i,:) = data_grp_psd.stats.std.pxx_period.loc{loc}.period{prd}.cond{cond};
                sum_array_grp_psd.pxx_f_period.mean(i,:) = data_grp_psd.stats.mean.pxx_f_period.loc{loc}.period{prd}.cond{cond};
                sum_array_grp_psd.pxx_f_period.std(i,:) = data_grp_psd.stats.std.pxx_f_period.loc{loc}.period{prd}.cond{cond};
                sum_array_grp_psd.pxx_period.col_freq{i,:} = round(sum_array_grp_psd.pxx_f_period.mean(i,:),1);                
                i = i +1;
            end
        end
    end
    
    %% STATS TEST - UNPAIRED T-TEST (is pxx(f) diff from zero?)
    % If psd is z-transformed and sham substracted it tests a difference from baseline that is beyond the difference of sham from baseline.
    p_threshold1 = 0.05;
    p_threshold2 = 0.05/2;
    i = 1;
    for loc=1:num_locations
        for prd=1:num_periods
            for cond=1:num_conditions-1   % 3rd cond is sham so ignored if sham substracted
                % stats test
                tmp = data_all_sbjcts_psd.pxx_period.loc{loc}.period{prd}.cond{cond};
                [h, p, ci, stats] = ttest(tmp);
                data_grp_psd.stats_test.ttest_p.pxx_period.loc{loc}.period{prd}(cond,:) = p;
                sig1 = p<p_threshold1;
                data_grp_psd.stats_test.ttest_p_sig1.pxx_period.loc{loc}.period{prd}(cond,:) = sig1;
                sig2 = p<p_threshold2;
                data_grp_psd.stats_test.ttest_p_sig2.pxx_period.loc{loc}.period{prd}(cond,:) = sig2;
                
                % put in a summary array
                sum_array_grp_psd.stats_test.ttest.row_dscrp{i,1} = [dscrp_loc{loc},'_', dscrp_prd{prd},'_',dscrp_cond{cond}];
                sum_array_grp_psd.stats_test.ttest.col_freq{i,:} = round(data_grp_psd.stats.std.pxx_period.loc{loc}.period{prd}.cond{cond},1);
                sum_array_grp_psd.stats_test.ttest.p(i,:) = p;
                sum_array_grp_psd.stats_test.ttest.sig1(i,:) = sig1;
                sum_array_grp_psd.stats_test.ttest.sig2(i,:) = sig2;
                
                i = i +1;
            end
        end
    end
    
    %% STATS TEST - PAIRED T-TEST (is pxx(f) of peak diff from trough?)
    p_threshold = 0.05;
    i = 1;
    for loc=1:num_locations
        for prd=2:num_periods
            % 3rd cond is sham so ignored if sham substracted
            tmp1 = data_all_sbjcts_psd.pxx_period.loc{loc}.period{prd}.cond{1}; % peak
            tmp2 = data_all_sbjcts_psd.pxx_period.loc{loc}.period{prd}.cond{2}; % trough
            [h, p, ci, stats] = ttest(tmp1, tmp2);
            data_grp_psd.stats_test.ttest2_p.pxx_period.loc{loc}.period{prd} = p;
            sig = p<p_threshold;
            data_grp_psd.stats_test.ttest2_p_sig.pxx_period.loc{loc}.period{prd} = sig;
            
            % put in a summary array
            sum_array_grp_psd.stats_test.ttest2.row_dscrp{i,1} = [dscrp_loc{loc},'_', dscrp_prd{prd}];
            sum_array_grp_psd.stats_test.ttest.col_freq{i,:} = round(data_grp_psd.stats.std.pxx_period.loc{loc}.period{prd}.cond{1},1);
            sum_array_grp_psd.stats_test.ttest2.p(i,:) = p;
            sum_array_grp_psd.stats_test.ttest2.sig(i,:) = sig;
            
            i = i +1;
        end
    end
    
    %% PLOT GROUP PSD
    if(plot_grp_psd)
        for loc=1:num_locations
            for prd=2:3 % EO and EC (basline is not plotted)
                h1 = figure;
                title = ['pwr_spctrm_grp_','_', dscrp_loc{loc},'_', dscrp_prd{prd}];
                set(h1,'name',title,'numbertitle','off');
                set(gca, 'FontSize', PlotFontSize, 'LineWidth',PlotAxesLineWidth); %<- Set properties
                for cond = 1:2
                    hold on
                    % extract vectors
                    v_f = data_grp_psd.stats.mean.pxx_f_period.loc{loc}.period{prd}.cond{cond};
                    v_mean = data_grp_psd.stats.mean.pxx_period.loc{loc}.period{prd}.cond{cond};
                    %                 v_err = data_grp_psd.stats.std.pxx_period.loc{loc}.period{prd}.cond{cond}; % std
                    v_err = data_grp_psd.stats.std.pxx_period.loc{loc}.period{prd}.cond{cond}/sqrt(data_grp_psd.stats.n.pxx_period.loc{loc}.period{prd}.cond{cond}); % sem
                    % plot mean
                    plot(v_f, v_mean ,'-o','color',ColorLine{cond},'LineWidth',2,'MarkerFaceColor',ColorLine{cond});
                    % plot std
                    errorbar(v_f, v_mean, v_err,'color',ColorLine{cond},'LineWidth',2);
                    hold off
                    xlim([0 25]);
                    max_y(cond) = max((v_mean)+(v_err))+0.01; max_y = max(max_y);
                    min_y(cond) = min((v_mean)-(v_err))-0.01; min_y = min(min_y);
                    ylim([min_y max_y]);
                    set(gca,'TickDir','out');
                    set(gca,'ticklength',3*get(gca,'ticklength'))
                    xlabel({'Frequency (Hz)'})
                    ylabel({'PSD (z-score)'})
                    if (saving_plot)
                        set(gcf,'PaperUnits','inches');
                        set(gcf,'PaperPosition',[0 0 7 6]);
                        set(gcf,'PaperPositionMode','manual')
                        %                         name = [title,'.fig'];
                        %                         saveas(h1,name);
                        name = [title];
                        print(name,'-dtiff',PlotPrintResolution);
                    end
                end
            end
        end
    end
    
    %% SAVE STRUCTURE
    if(save_grp_psd_structure)
        save(['data_grp_psd','.mat'],'data_grp_psd');
        save(['sum_array_grp_psd','.mat'],'sum_array_grp_psd');
    end
        
end


%%
cd(currentFolder)