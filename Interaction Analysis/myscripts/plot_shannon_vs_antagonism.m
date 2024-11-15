function fighandle = plot_shannon_vs_antagonism(subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix,fignum)


fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [100 100 700 350];


num_samples = size(subsampled_composition_matrix,1);
num_lineages = size(subsampled_composition_matrix,2);
num_subjects = numel(unique({subsampled_composition_table.Subject{:}}));
mean_antagonism_frequency = mean(subsampled_ZOI_matrix);

max_antag = zeros(num_samples,1);
shannon = zeros(num_samples,1);

for s=1:num_samples
    sample_idxs = find(subsampled_composition_matrix(s,:)>0);
    max_antag(s) = max(mean_antagonism_frequency(sample_idxs));
    for n=1:numel(sample_idxs)
        shannon(s) = shannon(s) - subsampled_composition_matrix(s,(sample_idxs(n)))*log(subsampled_composition_matrix(s,(sample_idxs(n))));
    end
end


subplot(1,2,1)
scatter(max_antag,shannon)
mdl = fitlm(max_antag,shannon);
hold on
plot(mdl)

[linear_rho,linear_p] = corr(max_antag,shannon,'Type','Pearson');
[spearman_rho,spearman_p] = corr(max_antag,shannon,'Type','Spearman');

xlabel('Max Inhibition frequency (vs. All Lineages)')
ylabel('Shannon Diversity')
title(['Correlation of Shannon diversity and antagonism' newline 'Linear rho = ' num2str(linear_rho) ' p = ' num2str(linear_p) newline 'Spearman rho = ' num2str(spearman_rho) ' p = ' num2str(spearman_p)])


subplot(1,2,2)
[sorted_antagonism,idxs] = sort(max_antag);
sorted_shannon = shannon(idxs);
bar(sorted_shannon);
hold on
bar(-sorted_antagonism);
xlabel('Sample')
ylabel('Frequency')
legend('Shannon','Max Anatag','Location','southwest')
xticks(1:num_samples)
xticklabels(subsampled_composition_table.SID(idxs))


