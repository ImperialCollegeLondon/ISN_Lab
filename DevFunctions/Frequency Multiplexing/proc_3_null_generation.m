clearvars
close all
clc

%% Load psd data
filename = 'eeg_1.mat';

load(filename);
load(['sleep_psd','_',char(regexp(filename,'[0-9]','match')),'.mat']);
load(['sleep_multiplex_',char(regexp(filename,'[0-9]','match')),'.mat']);

%% Single-Epoch Randomisation
f_width = 1;
N = 1419; % number of surrogate epochs to generate

tol = 0.03; % multiplexing tolerance

num_bins = floor(max(eeg_psd.freq)/f_width)+1;

surrogate_triplet_count = zeros(N,eeg_multiplex.nc);
surrogate_epoch = cell(N,eeg_multiplex.nc);


% find peaks distribution
pks_list = cell(1,eeg_multiplex.nc);
prob_dist = cell(1,eeg_multiplex.nc);
for ch = 1:eeg_multiplex.nc
    pks_list{ch} = cat(1, eeg_multiplex.pks_freq{ch,:});
    % fitting the pdf, width is chosen to best retain shape of the
    % histogram
    prob_dist{ch} = fitdist(pks_list{ch},'kernel','kernel','normal','width',0.05,'support',[0 30.01]);
end

% find number of peak distribution
npks_prob_dist = cell(1, eeg_multiplex.nc);
for ch = 1:eeg_multiplex.nc
    %     npks_prob_dist{ch} = fitdist(eeg_psd.npks(ch,:)','kernel','kernel','normal','width',0.1,'support',[min(eeg_psd.npks(ch,:))-1 max(eeg_psd.npks(ch,:))+1]);
    npks_prob_dist{ch} = fitdist(eeg_psd.npks(ch,:)','gamma');
end

df = round(1./(eeg.ep_sz(1)/eeg.Fs),1,'significant');
npks_range = [min(eeg_psd.npks,[],2) max(eeg_psd.npks,[],2)];
nc = eeg_psd.nc;

parfor i = 1:N
    if mod(i,N/10) == 0
        fprintf("Running iteration %d\n", i)
    end
    
    for ch = 1:nc
        
        % filling in the epoch
        npks = round(random(npks_prob_dist{ch}));
        npks = max(npks_range(ch,1),min(npks_range(ch,2),npks)); % caps the npks to the max and min from actual
        % discretise the peak positions
        generated_epoch = roundNearest(random(prob_dist{ch},npks,1), df);
        % make sure all peaks are unique
        generated_epoch = unique(generated_epoch);
        while length(generated_epoch) < npks
            n = npks - length(generated_epoch);
            tmp = roundNearest(random(prob_dist{ch},n,1), df);
            generated_epoch = unique([generated_epoch;tmp]);
        end
        
        % find multiplexing
        [triplet_count, ~, ~, ~] = multiplex_find(generated_epoch,generated_epoch, tol, max(eeg_psd.freq), 3);
        
        % get triplet count for epoch
        surrogate_triplet_count(i,ch) = sum(triplet_count);
        % save generated epochs
        surrogate_epoch{i,ch} = generated_epoch;
        
    end
end

%% Single-epoch Randomisation analysis

% generate mean triplet count per channel per iteration
surrogate_triplet_mean = mean(surrogate_triplet_count,1);
surrogate_triplet_std = std(surrogate_triplet_count,0,1);

actual_triplet = cellfun(@(x) sum(x(:,2)), eeg_multiplex.triplet_count,'un', false);
actual_triplet = cell2mat(actual_triplet);
actual_triplet_mean = mean(actual_triplet,2);
actual_triplet_std = std(actual_triplet,0,2);

% t-test with alpha = 0.01
for i = 1:eeg_multiplex.nc
    [ttest_result.H(i), ttest_result.p(i)] = ttest(actual_triplet(i,:),surrogate_triplet_mean(i),'Alpha',0.01,'Tail','right');
end

eeg_multiplex.ttest.mean = ttest_result;

% calculated total possible triplets for surrogate
[N, nc] = size(surrogate_epoch);
surrogate_possible_triplet = zeros(N,nc);
parfor i = 1:N
    for ch = 1:nc
        surrogate_possible_triplet(i,ch) = numPossibleTriplet(surrogate_epoch{i,ch}, 30);
    end
end

% calculated total possible triplet for actual
actual_possible_triplet = zeros(eeg_psd.nc, eeg_psd.nepc);
for ch = 1:eeg_psd.nc
    for epch = 1:eeg_psd.nepc
        actual_possible_triplet(ch,epch) = numPossibleTriplet(eeg_psd.pks_freq{ch,epch}, 30);
    end
end

% calculate percentage of triplet from total possible triplets
actual_triplet_percentage = actual_triplet./actual_possible_triplet * 100;
surrogate_triplet_percentage = surrogate_triplet_count./surrogate_possible_triplet * 100;

actual_triplet_percentage_mean = mean(actual_triplet_percentage,2);
actual_triplet_percentage_std = std(actual_triplet_percentage,0,2);
surrogate_triplet_percentage_mean = mean(surrogate_triplet_percentage);
surrogate_triplet_percentage_std = std(surrogate_triplet_percentage,0,1);

% t-test with alpha = 0.01
for i = 1:eeg_multiplex.nc
    [ttest_result.H(i), ttest_result.p(i)] = ttest(actual_triplet_percentage(i,:),mean(surrogate_triplet_percentage(:,i)),'Alpha',0.01,'Tail','right');
end

% ks-test with alpha = 0.01
for i = 1:eeg_multiplex.nc
    [ks_result.H(i), ks_result.p(i)] = kstest2(actual_triplet_percentage(i,:),surrogate_triplet_percentage(:,i),'Alpha',0.01,'Tail','larger');
end
    
eeg_multiplex.ttest.percentage = ttest_result;
eeg_multiplex.kstest = ks_result;

%% Single Epoch Mean Visualisation
hold on
X = categorical(string(eeg_multiplex.channels));
Y = [actual_triplet_mean surrogate_triplet_mean'];
colorA = [75,163,195]./256;
colorB = [23,86,118]./256;
b = bar(X,Y,'FaceColor','Flat');
b(1).FaceColor = colorA;
b(2).FaceColor = colorB;
ylabel('Mean number of triplets');
legend({'Actual','Surrogate'},'FontSize',12);
set(gca,'FontSize',15)

%% Single Epoch Percentage Visualisation
colorA = [75,163,195]./256;
colorB = [186,50,79]./256;

for ch = 1:eeg_multiplex.nc
    figure();
    hold on;
    histogram(surrogate_triplet_percentage(:,ch),'Normalization','pdf','FaceColor',colorA);
    histogram(actual_triplet_percentage(ch,:),'Normalization','pdf','FaceColor',colorB);
    xlabel('% of possible triplet','FontSize',15);
    ylabel('Probability','FontSize',15);
    legend({'Surrogate',string(eeg_multiplex.channels{ch})},'FontSize',12);
    set(gca,'FontSize',15)
    
    hold off    
end

figure()

hold on

X = categorical(string(eeg_multiplex.channels));
Y = [actual_triplet_percentage_mean surrogate_triplet_percentage_mean'];
b = bar(X,Y,'FaceColor','Flat');
b(1).FaceColor = colorB;
b(2).FaceColor = colorA;
ylabel('% of Possible Triplets');
legend({'Actual','Surrogate'},'FontSize',12);
set(gca,'FontSize',15)

hold off

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

parfor i = 1:N
    shuffled_pks_freq = eeg_multiplex.pks_freq(:,rand_seed(i,:));
    shuffled_multiplex = struct;
    
    shuffled_multiplex = duoEpochMultiplex(shuffled_multiplex,...
        eeg_multiplex.nc,eeg_multiplex.nepc,shuffled_pks_freq, max(eeg_psd.freq));
    
    shuffled_num_triplet = sumCellArray(shuffled_multiplex.duo_epoch.triplet_count);
    shuffled_mean_triplet(i,:) = mean(shuffled_num_triplet,2);
    
end

shuffled_std_triplet = std(shuffled_mean_triplet,1);

actual_num_triplet = sumCellArray(eeg_multiplex.duo_epoch.sum_count);
actual_mean_triplet = mean(actual_num_triplet,2);

%% peak pos pdf visualisation
ch = 4;
x = 0:0.05:max(pks_list{ch});
y = pdf(prob_dist{ch},x);
histogram(pks_list{ch},'Normalization','pdf','FaceColor',[0.2 0.5 0.7],'LineWidth',1);
hold on
plot(x,y,'r-','LineWidth',2);

xlabel('Peak frequency (Hz)','FontSize',15)
ylabel('Probability','FontSize',15)
legend({string(eeg_psd.channels{ch}),'Kernel fit'},'FontSize',12);
set(gca,'FontSize',14)
hold off

%% peak count pdf visualisation
ch = 4;
x = min(eeg_psd.npks(ch,:)):1:max(eeg_psd.npks(ch,:));
y = pdf(npks_prob_dist{ch},x);
histogram(eeg_psd.npks(ch,:),'Normalization','pdf','FaceColor',[208,221,215]./256);
hold on
plot(x,y,'-','Color',[89,78,54]./256,'LineWidth',2);
xlabel('No. of peaks','FontSize',15)
ylabel('Probability','FontSize',15)
legend({string(eeg_psd.channels{ch}),'Gamma fit'},'FontSize',12);
set(gca,'FontSize',14)
hold off

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