function [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix, lineages_represented] = subsample_composition(composition_table,interactions_by_lineage, comp_idxs)

% produces matched composition table, composition_matrix, and ZOI_matrix
% using only lineages present in the samples indicated by comp_idxs

subsampled_composition_table = composition_table(comp_idxs,:);
composition_matrix = subsampled_composition_table{:,4:end};
lineages_in_comp = find(sum(composition_matrix)>0); % find lineages present in composition

% subsample interaction structure
idxs = ismember([interactions_by_lineage.metadata.Lineage],lineages_in_comp);
interactions_by_lineage_subsample = subsample_interaction_structure(interactions_by_lineage,idxs);

% filter composition table
lineages_represented = [interactions_by_lineage_subsample.metadata.Lineage]; % find lineages represented in experiments
subsampled_composition_table = subsampled_composition_table(:,[1 2 3 lineages_represented+3]);
subsampled_composition_matrix = subsampled_composition_table{:,4:end};
subsampled_ZOI_matrix = interactions_by_lineage_subsample.ZOI_call; 
