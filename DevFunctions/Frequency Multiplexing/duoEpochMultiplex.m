function eeg_multiplex = duoEpochMultiplex(eeg_multiplex, nc,nepc,pks_freq,maxFreq, varargin)
%Duo Epoch Multiplex Summary of this function goes here
%   Detailed explanation goes here
eeg_multiplex.nepc = nepc;
eeg_multiplex.nc = nc;

eeg_multiplex.duo_epoch.triplet_count = cell(nc, nepc);
eeg_multiplex.duo_epoch.triplet_component = cell(nc, nepc);
eeg_multiplex.duo_epoch.harmonic_count = cell(nc, nepc);
eeg_multiplex.duo_epoch.harmonic_component = cell(nc, nepc);
eeg_multiplex.duo_epoch.generated_pks = cell(nc, nepc);
eeg_multiplex.duo_epoch.diff_count = cell(nc, nepc);
eeg_multiplex.duo_epoch.diff_component = cell(nc, nepc);

suppress = false;

if ~isempty(varargin)
    suppress = true;
end

for ch = 1:nc
    if ~suppress
        fprintf('Evaluating channel %d\n',ch);
    end
    
    for epch_2 = 2:nepc
        
        epch_1 = epch_2 - 1;
        
        tol = 0.03; % percentage tolerance
        
        [triplet_count, triplet_component, harmonic_count, harmonic_component, diff_count, diff_component] = multiplex_find(pks_freq{ch,epch_1},pks_freq{ch,epch_2}, tol, maxFreq, 3);
        
        % save to data structure
        eeg_multiplex.duo_epoch.triplet_count{ch,epch_2} = triplet_count;
        eeg_multiplex.duo_epoch.triplet_component{ch,epch_2} = triplet_component;
        eeg_multiplex.duo_epoch.harmonic_count{ch,epch_2} = harmonic_count;
        eeg_multiplex.duo_epoch.harmonic_component{ch,epch_2} = harmonic_component;
        eeg_multiplex.duo_epoch.diff_count{ch,epch_2} = diff_count;
        eeg_multiplex.duo_epoch.diff_component{ch,epch_2} = diff_component;
        
        % clear variables
        is_same_pks = [];
        same_pks = [];
        is_new_pks = [];
        generated_pks = [];
        
        % check if same peak occurs for both epoch with percentage tolerance
        [is_same_pks, same_pks] = ismembertol(pks_freq{ch,epch_2},pks_freq{ch,epch_1}, tol,'OutputAllIndices', 1);
        % set peaks with more than one match within tolerance to false. (assume
        % band splitting is new peaks)
        is_same_pks(cell2mat(cellfun(@(x) length(x) > 1, same_pks,'un',false))) = 0;
        is_new_pks = ~is_same_pks;
        
        % find peaks that are new and matches with peaks from multiplexing
        [~, generated_pks] = ismember(find(is_new_pks), find(triplet_count + diff_count + harmonic_count));
        %         [~, generated_pks] = ismember(find(is_new_pks), find(triplet_count(:) + harmonic_count(:)));
        generated_pks = generated_pks(generated_pks > 0);
        
        eeg_multiplex.duo_epoch.generated_pks{ch,epch_2} = generated_pks;
        eeg_multiplex.duo_epoch.is_new_pks{ch,epch_2} = is_new_pks;
        
    end
end
end

