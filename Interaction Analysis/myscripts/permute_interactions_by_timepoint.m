function [simulation_structure] = permute_interactions_by_timepoint(interactions_by_lineage,composition_table,shuffle_option,num_sims)

% groupby_option: metadata column to group by
% shuffle_option: rows, cols, or both to shuffle antagonizers and/or baits


formatSpec = '%.3f';
num_lineages = numel(interactions_by_lineage.metadata);
metadata_table = struct2table(interactions_by_lineage.metadata);
rep_lineages = [interactions_by_lineage.metadata.Lineage];
composition_matrix = composition_table{:,rep_lineages+3};
num_reps = size(composition_matrix,2);

% 
per_sample_same_group_interaction_freq = [];
per_sample_diff_group_interaction_freq = [];
same_group_interaction_vector = [];
diff_group_interaction_vector = [];
for s=1:27
    this_sample_same_group_interaction_vector = [];
    this_sample_diff_group_interaction_vector = [];
    in_sample = composition_matrix(s,:)>0;
    for m=1:num_reps
        for n=1:num_reps
            if m~=n
                if in_sample(m)&&in_sample(n)
                    this_sample_same_group_interaction_vector = [this_sample_same_group_interaction_vector; interactions_by_lineage.ZOI_call(m,n).*interactions_by_lineage.ZOI_AUC(m,n)];
                elseif in_sample(n) && ~in_sample(m)
                    this_sample_diff_group_interaction_vector = [this_sample_diff_group_interaction_vector; interactions_by_lineage.ZOI_call(m,n).*interactions_by_lineage.ZOI_AUC(m,n)];
                end
            end
        end
    end
    if isempty(this_sample_same_group_interaction_vector)
        this_sample_same_group_interaction_vector = [0];
    end


    per_sample_same_group_interaction_freq(s) = nnz(this_sample_same_group_interaction_vector)./numel(this_sample_same_group_interaction_vector);
    per_sample_diff_group_interaction_freq(s) = nnz(this_sample_diff_group_interaction_vector)./numel(this_sample_diff_group_interaction_vector);
    same_group_interaction_vector = [same_group_interaction_vector; this_sample_same_group_interaction_vector];
    diff_group_interaction_vector = [diff_group_interaction_vector; this_sample_diff_group_interaction_vector];
end

% Save observed data
simulation_structure.same_group_interaction_vector = same_group_interaction_vector;
simulation_structure.diff_group_interaction_vector = diff_group_interaction_vector;
simulation_structure.same_group_interaction_freq_sims = nan(num_sims,1);
simulation_structure.diff_group_interaction_freq_sims = nan(num_sims,1);
simulation_structure.per_sample_same_group_interaction_freq = per_sample_same_group_interaction_freq;
simulation_structure.per_sample_diff_group_interaction_freq = per_sample_diff_group_interaction_freq;


for x=1:num_sims
    switch shuffle_option
        case 'rows'
            idxs = randperm(num_reps);
            ZOI_call = interactions_by_lineage.ZOI_call(idxs, :);
            ZOI_AUC = interactions_by_lineage.ZOI_AUC(idxs, :);
        case 'cols'
            idxs = randperm(num_reps);
            ZOI_call = interactions_by_lineage.ZOI_call(:,idxs);
            ZOI_AUC = interactions_by_lineage.ZOI_AUC(:,idxs);

        case 'both'
            idxs = randperm(num_reps);
            ZOI_call = interactions_by_lineage.ZOI_call(idxs, :);
            ZOI_AUC = interactions_by_lineage.ZOI_AUC(idxs, :);
            idxs = randperm(num_reps);
            ZOI_call = ZOI_call(:,idxs);
            ZOI_AUC = ZOI_AUC(:,idxs);
        case 'random'
            idxs = reshape(randperm(num_reps*num_reps),num_reps,num_reps);
            ZOI_call = interactions_by_lineage.ZOI_call(idxs);
            ZOI_AUC = interactions_by_lineage.ZOI_AUC(idxs);
        otherwise
            ZOI_call = interactions_by_lineage.ZOI_call;
            ZOI_AUC = interactions_by_lineage.ZOI_AUC;
    end
    
    % iterate through each group
    same_group_interaction_vector = [];
    diff_group_interaction_vector = [];
    for s=1:27
        in_sample = composition_matrix(s,:)>0;
        for m=1:num_reps
            for n=1:num_reps
                if m~=n
                    if in_sample(m)&&in_sample(n)
                        same_group_interaction_vector = [same_group_interaction_vector; ZOI_call(m,n).*ZOI_AUC(m,n)];
                    elseif in_sample(n) && ~in_sample(m)
                        diff_group_interaction_vector = [diff_group_interaction_vector; ZOI_call(m,n).*ZOI_AUC(m,n)];
                    end
                end
            end
        end
        if isempty(same_group_interaction_vector)
            same_group_interaction_vector = [0];
        end
    end
    simulation_structure.same_group_interaction_freq_sims(x) = nnz(same_group_interaction_vector)./numel(same_group_interaction_vector);
    simulation_structure.diff_group_interaction_freq_sims(x) = nnz(diff_group_interaction_vector)./numel(diff_group_interaction_vector);
end

