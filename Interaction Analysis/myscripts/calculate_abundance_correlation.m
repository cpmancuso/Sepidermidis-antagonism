function [abundance_structure] = calculate_abundance_correlation(subsampled_composition_matrix,subsampled_ZOI_matrix)

% Note, this function assumes that composition matrix and ZOI matrix have
% been subsampled such that ZOI matrix only contains lineages that are
% present in composition matrix. 

% Note, compositions are not normalized, so composition may not add to 1
% due to lineages not represented in the experiment.


num_samples = size(subsampled_composition_matrix,1);
mean_antagonism_frequency = mean(subsampled_ZOI_matrix,1);
mean_sensitivity_frequency = mean(subsampled_ZOI_matrix,2)';
mean_abundance = sum(subsampled_composition_matrix)./sum(subsampled_composition_matrix>0);


abundance_vector = [];
antagonism_vector = [];
sensitivity_vector = [];
for s=1:num_samples
    lineages_in_sample = find(subsampled_composition_matrix(s,:)>0); %nonzero
    abundance_vector = [abundance_vector, subsampled_composition_matrix(s,lineages_in_sample)];
    antagonism_vector = [antagonism_vector, mean_antagonism_frequency(lineages_in_sample)];
    sensitivity_vector = [sensitivity_vector, mean_sensitivity_frequency(lineages_in_sample)];
end
abundance_vector = abundance_vector';
antagonism_vector = antagonism_vector';
sensitivity_vector = sensitivity_vector';

abundance_structure.abundance_vector = abundance_vector;
abundance_structure.antagonism_vector = antagonism_vector;
abundance_structure.sensitivity_vector = sensitivity_vector;
abundance_structure.mean_abundance = mean_abundance;
abundance_structure.mean_antagonism_frequency = mean_antagonism_frequency;
abundance_structure.mean_sensitivity_frequency = mean_sensitivity_frequency;

[abundance_structure.ant_rho_linear,abundance_structure.ant_p_linear] = corr(abundance_vector,antagonism_vector,'Type','Pearson');
[abundance_structure.ant_rho_spearman,abundance_structure.ant_p_spearman] = corr(abundance_vector,antagonism_vector,'Type','Spearman');
[abundance_structure.sen_rho_linear,abundance_structure.sen_p_linear] = corr(abundance_vector,sensitivity_vector,'Type','Pearson');
[abundance_structure.sen_rho_spearman,abundance_structure.sen_p_spearman] = corr(abundance_vector,sensitivity_vector,'Type','Spearman');


