function fighandle = plot_agr_and_phylo_composition(interactions_by_lineage_sepi,rep_lineages,composition_table,composition_matrix,fignum)

num_samples = size(composition_table,1);
num_lineages = size(composition_matrix,2);

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [100 100 700 350];

phylogroups = unique(cellstr({interactions_by_lineage_sepi.metadata.Phylogroup}));
agr_types = unique([interactions_by_lineage_sepi.metadata.Agr_Type]);


% phylogroups
phylo_comp = [];
for p=1:numel(phylogroups)
    idxs = find(contains([interactions_by_lineage_sepi.metadata.Phylogroup],phylogroups{p}));
    phylo_comp(:,p) = sum(composition_matrix(:,idxs),2);
end
% normalize
phylo_comp = phylo_comp./sum(phylo_comp,2);

subplot(1,2,1)
bar(phylo_comp,'Stacked')
title('Phylogroup composition')
legend(phylogroups,'location','northoutside','Orientation','horizontal')
ylim([0,1])
ylabel('Relative Abundance')
xticks([1:num_samples])
xticklabels(composition_table.SID)

% agr types
agr_comp = [];
for a=1:numel(agr_types)
    idxs = find([interactions_by_lineage_sepi.metadata.Agr_Type]==agr_types(a));
    agr_comp(:,a) = sum(composition_matrix(:,idxs),2);
end
% normalize
agr_comp = agr_comp./sum(agr_comp,2);

subplot(1,2,2)
bar(agr_comp,'Stacked')
title('AGR composition')
ylabel('Relative Abundance')

legend(num2str(agr_types'),'location','northoutside','Orientation','horizontal')
ylim([0,1])
xticks([1:num_samples])
xticklabels(composition_table.SID)