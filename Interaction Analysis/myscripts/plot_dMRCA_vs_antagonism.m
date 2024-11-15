% function fighandle = plot_dMRCA_vs_antagonism(interactions_by_lineage_sepi,dMRCA_filename,ant_option,fignum)
dMRCA_filename = 'dMRCA_table';
ant_option = 'sensitivity';

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [fighandle.Position(1) fighandle.Position(2) 700 350];

comp_idxs = 1:size(composition_table,1);
[subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix,lineages_represented] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
[abundance_structure] = calculate_abundance_correlation(subsampled_composition_matrix,subsampled_ZOI_matrix);

% Unpack structure
abundance_vector = abundance_structure.abundance_vector;
antagonism_vector = abundance_structure.antagonism_vector;
sensitivity_vector = abundance_structure.sensitivity_vector;
mean_abundance = abundance_structure.mean_abundance;
mean_antagonism_frequency = abundance_structure.mean_antagonism_frequency;
mean_sensitivity_frequency = abundance_structure.mean_sensitivity_frequency;

switch ant_option
    case 'antagonism'
        interaction_vector = antagonism_vector;
        mean_interaction_frequency = mean_antagonism_frequency;
    case 'sensitivity'
        interaction_vector = sensitivity_vector;
        mean_interaction_frequency = mean_sensitivity_frequency;
end

%% Load dMRCA table
molecular_clock = 4.5;
dMRCA_table = readtable(dMRCA_filename);
rowsToKeep = ismember(dMRCA_table.LineageNumber,lineages_represented);
dMRCA_table = dMRCA_table(rowsToKeep,:);
for m=1:size(dMRCA_table,1)
    dMRCA_vector = str2num(dMRCA_table.dMRCA_subject_ugenotypes{m});
    dMRCA_table.mean_dMRCA(m) = mean(dMRCA_vector);
    dMRCA_table.mean_tMRCA(m) = dMRCA_table.mean_dMRCA(m)./molecular_clock;
end

max_dMRCA_by_lineage = nan(size(mean_interaction_frequency));
median_dMRCA_by_lineage = nan(size(mean_interaction_frequency));
for n=1:numel(lineages_represented)
    lineage = lineages_represented(n);
    idxs = find(dMRCA_table.LineageNumber) == lineage;
    if ~isempty(idxs)
        max_dMRCA_by_lineage(n) = max(dMRCA_table.mean_dMRCA(idxs));
        median_dMRCA_by_lineage(n) = max(dMRCA_table.mean_dMRCA(idxs));
    end
end

%% Plot results

subplot(1,2,1)
scatter(mean_interaction_frequency,max_dMRCA_by_lineage)
mdl = fitlm(mean_interaction_frequency,max_dMRCA_by_lineage);
hold on
plot(mdl)

[linear_rho,linear_p] = corr(max_dMRCA_by_lineage',mean_interaction_frequency','Type','Pearson');
[spearman_rho,spearman_p] = corr(max_dMRCA_by_lineage',mean_interaction_frequency','Type','Spearman');

ylabel('Max dMRCA of Lineage')
switch ant_option
    case 'antagonism'
        xlabel('Antagonism frequency (vs. All Lineages)')
        title(['Correlation of lineage abundance and antagonism' newline 'Linear rho = ' num2str(linear_rho) ' p = ' num2str(linear_p) newline 'Spearman rho = ' num2str(spearman_rho) ' p = ' num2str(spearman_p)])
    case 'sensitivity'
        xlabel('Sensitivity frequency (vs. All Lineages)')
        title(['Correlation of lineage abundance and sensitivity' newline 'Linear rho = ' num2str(linear_rho) ' p = ' num2str(linear_p) newline 'Spearman rho = ' num2str(spearman_rho) ' p = ' num2str(spearman_p)])
end

subplot(1,2,2)
[sorted_interaction,idxs] = sort(mean_interaction_frequency);
sorted_dMRCA = max_dMRCA_by_lineage(idxs)./max(max_dMRCA_by_lineage(idxs));
bar(sorted_dMRCA);
hold on
bar(-sorted_interaction);
xlabel('Lineage')
ylabel('Frequency')
switch ant_option
    case 'antagonism'
        legend('Relative dMRCA','Antagonism','Location','southwest')
    case 'sensitivity'
        legend('Relative dMRCA','Sensitivity','Location','southwest')
end
xticks([])
set(gca,'box','off')
pbaspect([1 1 1])

    
