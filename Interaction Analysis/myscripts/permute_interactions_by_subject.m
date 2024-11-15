function [simulation_structure] = plot_per_sample_interaction_frequency_difference(interactions_by_lineage,composition_matrix,groupby_option,shuffle_option,num_sims,fignum)
fighandle = figure(fignum);
clf(fignum)
fighandle.Position=[fighandle.Position(1) fighandle.Position(2) 200 200]; %approx 2in sq

% groupby_option: metadata column to group by
% shuffle_option: rows, cols, or both to shuffle antagonizers and/or baits


formatSpec = '%.3f';
num_lineages = numel(interactions_by_lineage.metadata);
metadata_table = struct2table(interactions_by_lineage.metadata);


% 
for s=1:27
    rep_lineages = composition_matrix(s,:)>0;
    within_sample = 
