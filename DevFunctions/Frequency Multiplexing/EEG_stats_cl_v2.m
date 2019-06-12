%% npks stats

all_stats.npks_stats.epch.npks = eeg_multiplex.npks;

% Mean number of peaks across all epochs

for no_chnnls = 1:ch_end
    all_stats.npks_stats.channel.mean_npks(no_chnnls,1) = mean(eeg_multiplex.npks(no_chnnls,:));
end

% Max number of peaks across all epochs

for no_chnnls = 1:ch_end
    all_stats.npks_stats.channel.max_npks(no_chnnls,1) = max(eeg_multiplex.npks(no_chnnls,:));clc
end

% Min number of peaks across all epochs

for no_chnnls = 1:ch_end
    all_stats.npks_stats.channel.min_npks(no_chnnls,1) = min(eeg_multiplex.npks(no_chnnls,:));
end

% Number of triplets per epoch

triplet_count = zeros(size((eeg_multiplex.triplet_count),1),length(eeg_multiplex.triplet_count));
for no_chnnls = 1:ch_end
    for no_epchs = 1:epch_end
%         triplet_count(no_chnnls,no_epchs) = sum(eeg_multiplex.triplet_count{no_chnnls,no_epchs}(:,2)>0);       
        for no_pks = 1:length((eeg_multiplex.triplet_count{no_chnnls,no_epchs}(:,2)))
            if (eeg_multiplex.triplet_count{no_chnnls,no_epchs}(no_pks,2))>0
                all_stats.npks_stats.epch.triplet_count(no_chnnls,no_epchs) =  triplet_count(no_chnnls,no_epchs)+ 1;
            end
        end
    end
end

% Mean number of triplets across all epochs

for no_chnnls = 1:ch_end
    all_stats.npks_stats.channel.mean_triplet_count(no_chnnls,1) = mean(triplet_count(no_chnnls,:));
end

%% Percentage stats

% Triplets as a percentage of total number of possible triplets

for no_chnnls = 1:ch_end
    for no_epchs = 1:epch_end
        combi = combnk(1:(all_stats.npks_stats.epch.npks(no_chnnls, no_epchs)),3);
        all_stats.npks_stats.epch.no_possible_triplets(no_chnnls, no_epchs) = length(combi);
    end
end

% Mean possible number of triplets across all epochs

for no_chnnls = 1:ch_end
    all_stats.npks_stats.channel.mean_no_poss_triplets(no_chnnls,1) = mean(all_stats.npks_stats.epch.no_possible_triplets(no_chnnls, no_epchs));
end

for no_chnnls = 1:ch_end
    for no_epchs = 1:epch_end
        all_stats.percentage_stats.epch.percentage_triplets(no_chnnls, no_epchs) = ...
            (all_stats.npks_stats.epch.triplet_count(no_chnnls, no_epchs) / ...
            all_stats.npks_stats.epch.no_possible_triplets(no_chnnls, no_epchs) * 100);
    end
end

% Mean percentage and sd of peaks per epoch that are resulting from mixing across all epochs

for no_chnnls = 1:ch_end
    all_stats.percentage_stats.channel.triplet_percentage_stats.mean(no_chnnls,1) = mean(all_stats.percentage_stats.epch.percentage_triplets(no_chnnls,:));
    all_stats.percentage_stats.channel.triplet_percentage_stats.std(no_chnnls,1) = std(all_stats.percentage_stats.epch.percentage_triplets(no_chnnls,:));
end

%% Percentage of epochs with triplets

no_epchs_w_pks = 0;

for no_chnnls = 1:ch_end
    for no_epchs = 1:epch_end
        if (eeg_multiplex.npks(no_chnnls,no_epchs)>0)
            no_epchs_w_pks =  no_epchs_w_pks + 1;
        end
    end
end

no_epchs_w_triplet = 0;

for no_chnnls = 1:ch_end
    for no_epchs = 1:epch_end
        if (all_stats.npks_stats.epch.triplet_count(no_chnnls,no_epchs)>0)
            no_epchs_w_triplet =  no_epchs_w_triplet + 1;
        end
    end
end

percentage_epchs_w_triplet = no_epchs_w_triplet./no_epchs_w_pks * 100;

% Put data into structure 

all_stats.epchs_w_pks = no_epchs_w_pks;
all_stats.epchs_w_triplet = no_epchs_w_triplet;
all_stats.percentage_epchs_w_triplet = percentage_epchs_w_triplet;


