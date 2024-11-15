function [simulation_structure] = permute_interactions_by_group(interactions_by_lineage,groupby_option,shuffle_option,num_sims,fignum)
fighandle = figure(fignum);
clf(fignum)
fighandle.Position=[fighandle.Position(1) fighandle.Position(2) 200 200]; %approx 2in sq

% groupby_option: metadata column to group by
% shuffle_option: rows, cols, or both to shuffle antagonizers and/or baits


formatSpec = '%.3f';
num_lineages = numel(interactions_by_lineage.metadata);
metadata_table = struct2table(interactions_by_lineage.metadata);

%group by sample
[group, id] = findgroups(metadata_table.(groupby_option));

same_group_interaction_vector = [];
diff_group_interaction_vector = [];
num_reps = size(interactions_by_lineage.ZOI_call,1);

% iterate through each group
for m=1:num_reps
    for n=1:num_reps
        if m~=n
            if group(m)==group(n)
                same_group_interaction_vector = [same_group_interaction_vector; interactions_by_lineage.ZOI_call(m,n).*interactions_by_lineage.ZOI_AUC(m,n)];
            else
                diff_group_interaction_vector = [diff_group_interaction_vector; interactions_by_lineage.ZOI_call(m,n).*interactions_by_lineage.ZOI_AUC(m,n)];
            end
        end
    end
end
% Save observed data
simulation_structure.same_group_interaction_vector = same_group_interaction_vector;
simulation_structure.diff_group_interaction_vector = diff_group_interaction_vector;
simulation_structure.same_group_interaction_freq_sims = nan(num_sims,1);
simulation_structure.diff_group_interaction_freq_sims = nan(num_sims,1);


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
    for m=1:num_reps
        for n=1:num_reps
            if m~=n
                if group(m)==group(n)
                    same_group_interaction_vector = [same_group_interaction_vector; ZOI_call(m,n).*ZOI_AUC(m,n)];
                else
                    diff_group_interaction_vector = [diff_group_interaction_vector; ZOI_call(m,n).*ZOI_AUC(m,n)];
                end
            end
        end
    end
    simulation_structure.same_group_interaction_freq_sims(x) = nnz(same_group_interaction_vector)./numel(same_group_interaction_vector);
    simulation_structure.diff_group_interaction_freq_sims(x) = nnz(diff_group_interaction_vector)./numel(diff_group_interaction_vector);
end

