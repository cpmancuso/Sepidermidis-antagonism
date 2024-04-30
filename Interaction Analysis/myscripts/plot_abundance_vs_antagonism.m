function [fighandle] = plot_abundance_vs_antagonism(abundance_structure,ant_option,fignum)

% Note, this function assumes that composition matrix and ZOI matrix have
% been subsampled such that ZOI matrix only contains lineages that are
% present in composition matrix. 

% Note, compositions are not normalized, so composition may not add to 1
% due to lineages not represented in the experiment.


fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [fighandle.Position(1) fighandle.Position(2) 700 350];

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

subplot(1,2,1)
scatter(interaction_vector,abundance_vector)
mdl = fitlm(interaction_vector,abundance_vector);
hold on
plot(mdl)

[linear_rho,linear_p] = corr(abundance_vector,interaction_vector,'Type','Pearson');
[spearman_rho,spearman_p] = corr(abundance_vector,interaction_vector,'Type','Spearman');

ylabel('Abundance of Lineage (Each Sample)')
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
sorted_abundance = mean_abundance(idxs);
bar(sorted_abundance);
hold on
bar(-sorted_interaction);
xlabel('Lineage')
ylabel('Frequency')
switch ant_option
    case 'antagonism'
        legend('Relative Abundance','Antagonism','Location','southwest')
    case 'sensitivity'
        legend('Relative Abundance','Sensitivity','Location','southwest')
end
xticks([])
set(gca,'box','off')
pbaspect([1 1 1])
