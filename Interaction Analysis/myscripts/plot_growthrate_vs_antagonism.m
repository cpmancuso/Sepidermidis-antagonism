function [fighandle] = plot_growthrate_vs_antagonism(interactions_by_lineage_sepi,interactions_by_isolate,composition_table,growth_filename,fignum)

% Note, compositions are not normalized, so composition may not add to 1
% due to lineages not represented in the experiment.

% Options
plot_option = 'median'; % 'median' or 'allpoints'

comp_idxs = 1:size(composition_table,1);
[subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix,lineages_represented] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
[abundance_structure] = calculate_abundance_correlation(subsampled_composition_matrix,subsampled_ZOI_matrix);

% Unpack structure
mean_antagonism_frequency = abundance_structure.mean_antagonism_frequency;
mean_sensitivity_frequency = abundance_structure.mean_sensitivity_frequency;

growth_table = readtable(growth_filename);

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [100 100 700 350];

% Parse growth rate for each lineage in experiment
num_lineages = numel(lineages_represented);
%set up default points for legend
switch plot_option
    case 'allpoints'
        scatter([-1],[1],10,'k','Filled','MarkerEdgeAlpha',0.5)
        hold on
        scatter([-1],[1],20,'k','Filled','MarkerEdgeAlpha',0.5)
        scatter([-1],[1],30,'k','Filled','MarkerEdgeAlpha',0.5)
        legend('Day 1','Day 2','Day 3','Location','north','Orientation','horizontal')
end
lineage_median_growth = [];
lineage_std_growth = [];
for n=1:num_lineages
    hold on
    lineage = lineages_represented(n);
    idxs = ([interactions_by_isolate.metadata.Lineage] == lineage);
    isolate_list = unique([interactions_by_isolate.metadata(idxs).Stock]);
    replicate_medians = [];
    replicate_stds = [];
    for i=1:numel(isolate_list)
        isolate = isolate_list(i);
        [replicate_medians(i),replicate_stds(i),growthrates,days] = get_lineage_growthrate(growth_table,isolate);
        switch plot_option
            case 'allpoints'
                scatter(n*ones(size(growthrates)),growthrates,days*20,'Filled','MarkerFaceAlpha',0.5,'HandleVisibility','off')
            case 'median'
                jitter = 0.25*(rand-0.5);
                % scatter(n,median(growthrates),10,'Filled','MarkerFaceAlpha',0.8,'HandleVisibility','off')
                errorbar(n+jitter,median(growthrates),std(growthrates)/sqrt(length(growthrates)),'o','MarkerFaceColor','auto','MarkerEdgeColor','auto','MarkerSize',5,'CapSize',0,'HandleVisibility','off')
        end
    end
    lineage_median_growth(n) = nanmedian(replicate_medians);
    lineage_std_growth(n) = sqrt(sum(replicate_stds.^2)) / numel(replicate_stds);
end
xlim([0 num_lineages+1])
ylim([0 1.2])
xticks(1:num_lineages)
xticklabels(lineages_represented)
xlabel('Lineage')
ylabel('Growth Rate [h-1]')

% plot antagonism vs growth rate
fighandle = figure(fignum+1);
clf(fignum+1)
fighandle.Position = [100 100 700 350];
median_interaction_frequency = mean_antagonism_frequency;
subplot(1,2,1)
scatter(median_interaction_frequency,lineage_median_growth)
mdl = fitlm(median_interaction_frequency,lineage_median_growth);
hold on
plot(mdl)

plot_idxs =  ~isnan(lineage_median_growth) & ~isnan(median_interaction_frequency);
median_growth_to_plot = lineage_median_growth(plot_idxs)';
median_interaction_frequency_to_plot = median_interaction_frequency(plot_idxs)';

[linear_rho,linear_p] = corr(median_growth_to_plot,median_interaction_frequency_to_plot,'Type','Pearson');
[spearman_rho,spearman_p] = corr(median_growth_to_plot,median_interaction_frequency_to_plot,'Type','Spearman');

ylabel('Lineage Growth Rate')
xlabel('Antagonism frequency (vs. All Lineages)')
title(['Correlation of growth rate and antagonism' newline 'Linear rho = ' num2str(linear_rho) ' p = ' num2str(linear_p) newline 'Spearman rho = ' num2str(spearman_rho) ' p = ' num2str(spearman_p)])


subplot(1,2,2)
[sorted_interaction,idxs] = sort(median_interaction_frequency_to_plot);
sorted_growth = median_growth_to_plot(idxs);
bar(sorted_growth);
hold on
bar(-sorted_interaction);
xlabel('Lineage')
ylabel('Frequency')
legend('Growth Rate','Antagonism','Location','southwest')


% plot sensitivity vs growth rate
fighandle = figure(fignum+2);
clf(fignum+2)
fighandle.Position = [100 100 700 350];
median_interaction_frequency = mean_sensitivity_frequency;
subplot(1,2,1)
scatter(median_interaction_frequency,lineage_median_growth)
mdl = fitlm(median_interaction_frequency,lineage_median_growth);
hold on
plot(mdl)

plot_idxs =  ~isnan(lineage_median_growth) & ~isnan(median_interaction_frequency);
median_growth_to_plot = lineage_median_growth(plot_idxs)';
median_interaction_frequency_to_plot = median_interaction_frequency(plot_idxs)';

[linear_rho,linear_p] = corr(median_growth_to_plot,median_interaction_frequency_to_plot,'Type','Pearson');
[spearman_rho,spearman_p] = corr(median_growth_to_plot,median_interaction_frequency_to_plot,'Type','Spearman');

ylabel('Lineage Growth Rate')
xlabel('Sensitivity frequency (vs. All Lineages)')
title(['Correlation of growth rate and sensitivity' newline 'Linear rho = ' num2str(linear_rho) ' p = ' num2str(linear_p) newline 'Spearman rho = ' num2str(spearman_rho) ' p = ' num2str(spearman_p)])


subplot(1,2,2)
[sorted_interaction,idxs] = sort(median_interaction_frequency_to_plot);
sorted_growth = median_growth_to_plot(idxs);
bar(sorted_growth);
hold on
bar(-sorted_interaction);
xlabel('Lineage')
ylabel('Frequency')
legend('Growth Rate','Sensitivity','Location','southwest')


%% Subfunctions
    function [median_growthrate,std_growthrate,growthrates,days] = get_lineage_growthrate(growth_table,isolate_list)
        R_threshold = 0.99;
        [growth_idxs] = ismember(growth_table.Culture,isolate_list);
        growth_idxs = logical(growth_idxs.*(growth_table.Rsq>R_threshold));
        
        if sum(growth_idxs)
            growthrates = growth_table.GrowthRate(growth_idxs);
            days = growth_table.Day(growth_idxs);
            median_growthrate = median(growthrates);
            std_growthrate = median(growthrates);
        else
            growthrates = nan;
            days = nan;
            median_growthrate = nan;
            std_growthrate = nan;
        end
    end
end