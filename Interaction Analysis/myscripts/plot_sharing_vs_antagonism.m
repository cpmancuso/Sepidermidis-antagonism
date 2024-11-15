function fighandle = plot_sharing_vs_antagonism(subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix,fignum)

% Note, this function assumes that composition matrix and ZOI matrix have
% been subsampled such that ZOI matrix only contains lineages that are
% present in composition matrix. 

% Note, compositions are not normalized, so composition may not add to 1
% due to lineages not represented in the experiment.

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [100 100 700 350];


num_samples = size(subsampled_composition_matrix,1);
num_lineages = size(subsampled_composition_matrix,2);
num_subjects = numel(unique({subsampled_composition_table.Subject{:}}));
mean_antagonism_frequency = mean(subsampled_ZOI_matrix);

for n=1:num_lineages
    subject_idxs = find(subsampled_composition_matrix(:,n)>0);
    num_subjects_colonized(n) = numel(unique({subsampled_composition_table.Subject{subject_idxs}}));
end
mean_antagonism_frequency = mean_antagonism_frequency';
num_subjects_colonized = num_subjects_colonized'./num_subjects;

subplot(1,2,1)
scatter(mean_antagonism_frequency,num_subjects_colonized)
mdl = fitlm(mean_antagonism_frequency,num_subjects_colonized);
hold on
plot(mdl)

[linear_rho,linear_p] = corr(mean_antagonism_frequency,num_subjects_colonized,'Type','Pearson');
[spearman_rho,spearman_p] = corr(mean_antagonism_frequency,num_subjects_colonized,'Type','Spearman');

xlabel('Inhibition frequency (vs. All Lineages)')
ylabel('Fraction of subjects with lineage')
title(['Correlation of lineage sharing and antagonism' newline 'Linear rho = ' num2str(linear_rho) ' p = ' num2str(linear_p) newline 'Spearman rho = ' num2str(spearman_rho) ' p = ' num2str(spearman_p)])


subplot(1,2,2)
[sorted_antagonism,idxs] = sort(mean_antagonism_frequency);
sorted_num_subjects_colonized = num_subjects_colonized(idxs);
bar(sorted_num_subjects_colonized);
hold on
bar(-sorted_antagonism);
xlabel('Lineage')
ylabel('Frequency')
legend('Fraction of Lineages','Antagonism','Location','southwest')
ylim([-1 1])



