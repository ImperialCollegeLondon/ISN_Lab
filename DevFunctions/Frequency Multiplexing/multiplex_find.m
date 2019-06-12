function [triplet_num,triplet_component, harmonic_num, harmonic_component, varargout] = multiplex_find(epch_1,epch_2, tol, maxF, max_harmonic)
%MULTIPLEX_FIND From two vectors of frequencies, count the number of
%triplets and multiplexing peaks.
%   epch_1 : first epoch
%   epch_2 : second epoch
%   tol    : fractional tolerance
%   maxF   : maximum frequency
%   max_harmonic : maximum order of harmonics considered
%
%   triplet_num        : number of triplets contributing
%   triplet_component  : frequency of the two other peaks in the triplet
%   harmonic_num       : number of peaks' harmonic contributing
%   harmonic_component : the fundamental frequency of the peak
%
%   OPTIONAL
%   diff_count : number of peak differences contributing
%   diff_component : frequency of the parent peaks

nout = max(nargout,0);
find_diff = nout > 0;

npks_1 = length(epch_1);
npks_2 = length(epch_2);

% Clear variables
if find_diff
    pks_freq_diff = zeros(npks_1);
    pks_freq_closest_pks_freq_diff = zeros(npks_1);
    pks_freq_closest_pks_freq_diff_indx = zeros(npks_1);
    pks_freq_diff_exist = zeros(npks_1);
    diff_count = zeros(1,npks_2);
    diff_component = {};
end
pks_freq_sum = zeros(npks_1);
pks_freq_closest_pks_freq_sum = zeros(npks_1);
pks_freq_closest_pks_freq_sum_indx = zeros(npks_1);
pks_freq_sum_exist = zeros(npks_1);
pks_freq_harmonic = zeros(npks_1, max_harmonic);
pks_freq_closest_pks_freq_harmonic = pks_freq_harmonic;
pks_freq_closest_pks_freq_harmonic_indx = pks_freq_harmonic;
pks_freq_harmonic_exist = pks_freq_harmonic;

triplet_num = zeros(1,npks_2);
triplet_component = {};
harmonic_num = zeros(1,npks_2);
harmonic_component = {};

% Multiplexing and Harmonic analysis
for i = 1:npks_1
    % multiplexing
    for k = i+1:npks_1
        % diff
        % calculate frequency difference
        if find_diff
            pks_freq_diff(i,k) = abs(epch_1(k) - epch_1(i));
            % find index of the peak that matches with the frequency difference
            pks_freq_closest_pks_freq_diff_indx(i,k) = closest(pks_freq_diff(i,k), epch_2);
            % use the index to find the actual peak frequency
            pks_freq_closest_pks_freq_diff(i,k) = epch_2(pks_freq_closest_pks_freq_diff_indx(i,k));
            % check if the actual peak and the calculated peak are within tolerance
            pks_freq_diff_exist(i,k) = (abs(pks_freq_closest_pks_freq_diff(i,k) - pks_freq_diff(i,k)) < tol * pks_freq_diff(i,k));
            % reject peak matches that are its parent peaks
            if(pks_freq_closest_pks_freq_diff_indx(i,k) == i || pks_freq_closest_pks_freq_diff_indx(i,k) == k)
                pks_freq_diff_exist(i,k) = 0;
                pks_freq_closest_pks_freq_diff_indx(i,k) = 0;
                pks_freq_closest_pks_freq_diff(i,k) = 0;
            end
        end
        
        % sum
        % calculate frequency sum
        pks_freq_sum(i,k) = abs(epch_1(k) + epch_1(i));
        % find index of the peak that matches with the frequency sum
        pks_freq_closest_pks_freq_sum_indx(i,k) = closest(pks_freq_sum(i,k), epch_2);
        % use the index to find the actual peak frequency
        pks_freq_closest_pks_freq_sum(i,k) = epch_2(pks_freq_closest_pks_freq_sum_indx(i,k));
        % check if the actual peak and the calculated peak are within tolerance
        pks_freq_sum_exist(i,k) = (abs(pks_freq_closest_pks_freq_sum(i,k) - pks_freq_sum(i,k)) < tol * pks_freq_sum(i,k));
        % reject peak matches that are its parent peaks
        if(pks_freq_closest_pks_freq_sum_indx(i,k) == i || pks_freq_closest_pks_freq_sum_indx(i,k) == k)
            pks_freq_sum_exist(i,k) = 0;
            pks_freq_closest_pks_freq_sum_indx(i,k) = 0;
            pks_freq_closest_pks_freq_sum(i,k) = 0;
        end
    end
    
    % harmonics
    % calculate maximum multiple before exceeding frequency range
    mul = floor(maxF/epch_1(i));
    if mul >= 2
        for k = 2:max_harmonic
            % calculate the k-th harmonic
            pks_freq_harmonic(i,k) = epch_1(i) * k;
            % find index of the peak that matches with the harmonic
            pks_freq_closest_pks_freq_harmonic_indx(i,k) = closest(pks_freq_harmonic(i,k), epch_2);
            % use the index to find the actual peak frequency
            pks_freq_closest_pks_freq_harmonic(i,k) = epch_2(pks_freq_closest_pks_freq_harmonic_indx(i,k));
            % check if the actual peak and the calculated peak are within tolerance
            pks_freq_harmonic_exist(i,k) = (abs(pks_freq_closest_pks_freq_harmonic(i,k) - pks_freq_harmonic(i,k)) < tol * pks_freq_harmonic(i,k));
            % reject peak matches that are its fundamental frequency
            if(pks_freq_closest_pks_freq_harmonic_indx(i,k) == i)
                pks_freq_harmonic_exist(i,k) = 0;
                pks_freq_closest_pks_freq_harmonic_indx(i,k) = 0;
                pks_freq_closest_pks_freq_harmonic(i,k) = 0;
            end
        end
    else
        pks_freq_harmonic(i,1) = 0;
    end
end

% Reset elements that correspond to non-existent
% multiplexing/harmonics
pks_freq_closest_pks_freq_sum_indx(~logical(pks_freq_sum_exist)) = 0;
pks_freq_closest_pks_freq_sum(~logical(pks_freq_sum_exist)) = 0;
% pks_freq_closest_pks_freq_diff(~logical(pks_freq_diff_exist)) = 0;
% pks_freq_closest_pks_freq_diff_indx(~logical(pks_freq_diff_exist)) = 0;
pks_freq_closest_pks_freq_harmonic(~logical(pks_freq_harmonic_exist)) = 0;
pks_freq_closest_pks_freq_harmonic_indx(~logical(pks_freq_harmonic_exist)) = 0;

% Summarise mulitplexing contributions
for i = 1:npks_2
    % count sum
    triplet_num(i) = sum(sum(pks_freq_closest_pks_freq_sum_indx == i));
    if triplet_num(i) > 0
        % find index of parent peaks
        tmp = find(pks_freq_closest_pks_freq_sum_indx == i);
        [x,y] = ind2sub(size(pks_freq_closest_pks_freq_sum_indx), tmp);
        for j = 1:triplet_num(i)
            triplet_component{i,j} = [x(j) y(j)];
        end
    else
        triplet_component{i,1} = [];
    end
    
    % count diff
    if find_diff
        diff_count(i) = sum(sum(pks_freq_closest_pks_freq_diff_indx == i));
        if diff_count(i) > 0
            % find index of parent peaks
            tmp = find(pks_freq_closest_pks_freq_diff_indx == i);
            [x,y] = ind2sub(size(pks_freq_closest_pks_freq_diff_indx), tmp);
            for j = 1:diff_count(i)
                diff_component{i,j} = [x(j) y(j)];
            end
        else
            diff_component{i,1} = [];
        end
    end
    
    % count harmonics
    harmonic_num(i) = sum(sum(pks_freq_closest_pks_freq_harmonic_indx == i));
    if harmonic_num(i) > 0
        % find index of parent peak and order of harmonic
        tmp = find(pks_freq_closest_pks_freq_harmonic_indx == i);
        [x,y] = ind2sub(size(pks_freq_closest_pks_freq_harmonic_indx), tmp);
        for j = 1:harmonic_num(i)
            harmonic_component{i,j} = [x(j) y(j)];
        end
    else
        harmonic_component{i,1} = [];
    end
end

if find_diff
    varargout(1) = {diff_count};
    varargout(2) = {diff_component};
end

end