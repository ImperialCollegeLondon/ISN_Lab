function [cond_order, rand_cond_order] = ISN_ConditionBalancing(num_cond)
%ISN_ConditionBalancing Creates a condition order and randomised condition
%order for any number of conditions (num_cond)
%   [cond_order, rand_cond_order] = ISN_ConditionBalancing(num_cond)
%
% Eddy Rhodes, July 2019
%
% Ref: Sharma, V.K. (1975). An Easy Method of Constructing Latin Square 
% Designs Balanced for the Immediate Residualand Other Order Effects.
% Canadian Journal of Statistics. 3(1). p 119-124
%% DETERMINE IF THE NUMBER OF CONDITIONS IS EVEN OR ODD
% This is required to allow for the balanced square (squares if odd) to be
% correctly filled
odd_n = mod(num_cond,2);
if odd_n ==1
    odd_i = 1:2:num_cond;
    even_i = 2:2:num_cond-1;
    odd_l = length(odd_i);
    even_l = length(even_i);
else
    odd_i = 1:2:num_cond-1;
    even_i = 2:2:num_cond;
    odd_l = length(odd_i);
    even_l = length(even_i);
end
%% CREATE 1st ROW OF BALANCED SQUARE
row1(odd_i) = 1:odd_l;
row1(even_i) = num_cond:-1:num_cond-(even_l-1);
if odd_n == 1
    row2(even_i) = 1:even_l;
    row2(odd_i) = num_cond:-1:num_cond-(odd_l-1);
end

%% FILL REST OF THE BALANCED SQUARE(S)
sq1 = zeros(num_cond,num_cond);
sq1(1,:) = row1;
for n = 2:num_cond
    for i = 1:num_cond
        sq1(n,i) = sq1(n-1,i)+1;
        if sq1(n,i) > num_cond
            sq1(n,i) = 1;
        end
    end
end
if odd_n == 1
    sq2 = zeros(num_cond,num_cond);
    sq2(1,:) = row2;
    for n = 2:num_cond
        for i = 1:num_cond
            sq2(n,i) = sq2(n-1,i)+1;
            if sq2(n,i) > num_cond
                sq2(n,i) = 1;
            end
        end
    end
end
tmp = sq1;
if odd_n == 1
    tmp = [tmp;sq2];
end
%% CHECK THAT EACH COLUMN IS COUNTERBALANCED
h1 = figure;
figtitle = 'Check_counterbalance_in_each_column';
set(h1,'name',figtitle,'numbertitle','off');
y = floor(sqrt(num_cond));
x = ceil(num_cond/y);
for n = num_cond:-1:1
    subplot(y,x,n)
    histogram(tmp(:,n))
    xlabel('Stim Condition')
    ylabel('No. of Repetitions')
    title(['Column no. ' num2str(n)])
    % STORE CHECK VALUE
    tmp_u = unique(tmp(:,n));
    for u = size(tmp_u,1):-1:1
        tmp_ui(u) = size(find(tmp(:,n)==tmp_u(u)),1);
    end
    tmp_chk(n) = sum(diff(tmp_ui)~=0);
end
if sum(tmp_chk)>0
    warning('Each column is not counterbalanced')
end
clear tmp_chk tmp_ui tmp_u
%% CHECK THAT THE ORDER IS COUNTERBALANCED
tmp_c = cell(num_cond,1);
for n = 1:(size(tmp,2)-1)
    for m = 1:size(tmp,1)
        tmpp = tmp(m,n);
        tmppp = tmp(m,n+1);
        tmp_c{tmpp} = [tmp_c{tmpp},tmppp];
    end
end
h1 = figure;
figtitle = 'Check_cond_order_is_not_biased';
set(h1,'name',figtitle,'numbertitle','off');
for n = 1:num_cond
    subplot(y,x,n)
    histogram(tmp_c{n})
    xlabel('Stim Condition')
    ylabel('No. of Repetitions')
    title(['Condition no. ' num2str(n)])
    % STORE CHECK VALUE
    tmp_u = unique(tmp_c{n});
    for u = size(tmp_u,2):-1:1
        tmp_ui(u) = size(find(tmp_c{n}==tmp_u(u)),2);
    end
    tmp_chk(n) = sum(diff(tmp_ui)~=0);
end
if sum(tmp_chk)>0
    warning('There is an order bias')
end
%%
cond_order = tmp;
ind = 1:size(cond_order,1);
rand = randperm(size(cond_order,1),size(cond_order,1));
ind=ind(rand);
rand_cond_order = cond_order(rand,:);
end

