clearvars
close all
clc

%% Load psd data
filename = 'eeg_1.mat';

load(['sleep_psd','_',char(regexp(filename,'[0-9]','match')),'.mat']);
load(['sleep_multiplex_',char(regexp(filename,'[0-9]','match')),'.mat']);

%% Single-Epoch Randomisation
f_width = 1;
N = 1000; % number of randomisation

tol = 0.03;

num_bins = floor(max(eeg_psd.freq)/f_width)+1;

shuffle_triplet_count = zeros(eeg_multiplex.nc,eeg_multiplex.nepc);
shuffle_triplet_mean = zeros(eeg_multiplex.nc,N);

% cells of peak positions
pks_list = cell(eeg_psd.nc,num_bins);
pks_edge = [0:f_width:max(eeg_psd.freq)+1];
pks_count = cell(eeg_psd.nc,eeg_psd.nepc);

% if unique peak position is used for each epoch
pks_unique = true;
if pks_unique
    pks_used = cell(eeg_psd.nc,num_bins);
end

% collect list of peaks positions for each bin
for ch = 1:eeg_multiplex.nc
    for epch = 1:eeg_multiplex.nepc
        for pk = 1:num_bins
            pks_count{ch,epch}(pk) = sum(pks_edge(pk) <= eeg_multiplex.pks_freq{ch,epch} & ...
                eeg_multiplex.pks_freq{ch,epch} < pks_edge(pk+1));
            
            pks_list{ch,pk} =[pks_list{ch,pk}; eeg_multiplex.pks_freq{ch,epch}(...
                pks_edge(pk) <= eeg_multiplex.pks_freq{ch,epch} & ...
                eeg_multiplex.pks_freq{ch,epch} < pks_edge(pk+1))];
            
            
            if pks_unique
                pks_used{ch,pk} = zeros(1,length(pks_list{ch,pk}));
            end
        end
    end
end

pks_list_mean = cellfun(@mean, pks_list);
pks_list_std = cellfun(@std, pks_list);

for i = 1:N
    if mod(i,10) == 0
        fprintf("Running iteration %d\n", i)
    end
    
    pks_used = cellfun(@(x) zeros(1,length(x)), pks_used, 'un', false);
    generated_epochs = cell(eeg_psd.nc,eeg_psd.nepc);
    for ch = 1:eeg_multiplex.nc
        for epch = 1:eeg_multiplex.nepc
            generated_epochs{ch,epch} = zeros(1,eeg_psd.npks(ch,epch));
            
            pk_idx = 1;
            
            % generating random epochs
            for pk = 1:num_bins
                if pks_count{ch,epch}(pk)
                    for j = 1:pks_count{ch,epch}(pk)
                        
                        rand_pk = randi(length(pks_list{ch,pk}),1);
                        
                        if pks_unique
                            while pks_used{ch,pk}(rand_pk)
                                rand_pk = randi(length(pks_list{ch,pk}),1);
                            end
                        end
                        
                        generated_epochs{ch,epch}(pk_idx) = pks_list{ch,pk}(rand_pk);
                        pks_used{ch,pk}(rand_pk) = true;
                        pk_idx = pk_idx + 1;
                    end
                end
            end
            
            % find multiplexing
            [triplet_count, ~, ~, ~] = multiplex_find(generated_epochs{ch,epch},generated_epochs{ch,epch}, tol, max(eeg_psd.freq), 3);
            
            % get triplet count for epoch            
            shuffle_triplet_count(ch,epch) = sum(triplet_count);
                      
        end
    end
    
    % generate mean triplet count per channel per iteration
    shuffle_triplet_mean(:,i) = mean(shuffle_triplet_count,2);  
end

%% Single-epoch Randmoisation analysis

shuffle_triplet_mean_mean = mean(shuffle_triplet_mean,2);
shuffle_triplet_mean_std = std(shuffle_triplet_mean,0,2);

actual_triplet = cellfun(@(x) sum(x(:,2)), eeg_multiplex.triplet_count,'un', false);
actual_triplet = cell2mat(actual_triplet);
actual_triplet_mean = mean(actual_triplet,2);
actual_triplet_std = std(actual_triplet,0,2);

% t-test with alpha = 0.01
for i = 1:eeg_multiplex.nc
    [ttest_result.H(i), ttest_result.p(i)] = ttest(actual_triplet(i,:),shuffle_triplet_mean_mean(i),'Alpha',0.01,'Tail','right');
end

eeg_multiplex.ttest = ttest_result;

%% Duo-epoch Randomisation
N = 1000;
shuffled_mean_triplet = zeros(N,eeg_multiplex.nc);

rand_seed = zeros(N,eeg_multiplex.nepc);
rand_seed(1,:) = randperm(eeg_multiplex.nepc);

for i = 2:N
    tmp = randperm(eeg_multiplex.nepc);
    while duplicatePerm(rand_seed(1:i-1,:), tmp)
        tmp = randperm(eeg_multiplex.nepc);
    end
    rand_seed(i,:) = tmp;
end

for i = 1:N
    shuffled_pks_freq = eeg_multiplex.pks_freq(:,rand_seed(i,:));
    shuffled_multiplex = struct;
    
    shuffled_multiplex = duoEpochMultiplex(shuffled_multiplex,...
        eeg_multiplex.nc,eeg_multiplex.nepc,shuffled_pks_freq, max(eeg_psd.freq));
    
    shuffled_num_triplet = sumCellArray(shuffled_multiplex.duo_epoch.triplet_count);
    shuffled_mean_triplet(i,:) = mean(shuffled_num_triplet,2);
    
end

actual_possible_tri_gen_pairs = zeros(eeg_multiplex.nc,eeg_multiplex.nepc);

for ch = 1:eeg_multiplex.nc
    for epch = 1:eeg_multiplex.nepc
        combi = combnk(eeg_multiplex.pks_freq{ch,epch},2); % number of pairs
        
        actual_possible_tri_gen_pairs(ch, epch) = length(combi);
    end
end

shuffled_std_triplet = std(shuffled_mean_triplet,1);

actual_num_triplet = sumCellArray(eeg_multiplex.duo_epoch.sum_count);
actual_mean_triplet = mean(actual_num_triplet,2);

actual_percent_triplet = actual_num_triplet(:,2:end)./actual_possible_tri_gen_pairs(:,1:end-1);
actual_percent_triplet = [zeros(4,1) actual_percent_triplet];

%% Utilities functions
function summed_output = sumCellArray(cArray)

cArray_sz = size(cArray);
summed_output = zeros(cArray_sz);

for ch = 1:cArray_sz(1)
    for epch = 1:cArray_sz(2)
        summed_output(ch,epch) = sum(cArray{ch,epch});
    end
end

end

function output = duplicatePerm(A,B)

[row,~] = size(A);

for i = 1:row
    if isempty(A(i,:) -B)
        output = 1;
        return
    end
end

output = 0;

end