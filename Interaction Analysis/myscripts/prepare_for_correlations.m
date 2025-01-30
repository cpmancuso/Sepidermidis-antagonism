function correlation_structure = prepare_for_correlations(subsampled_composition_table,subsampled_ZOI_matrix,growth_filename,interactions_by_isolate,lineage_labels);

% Note, this function assumes that composition matrix and ZOI matrix have
% been subsampled such that ZOI matrix only contains lineages that are
% present in composition matrix. 

% Note, compositions are not normalized, so composition may not add to 1
% due to lineages not represented in the experiment.

% Load composition matrix and growth table
subsampled_composition_matrix = subsampled_composition_table{:,4:end};
growth_table = readtable(growth_filename);

num_samples = size(subsampled_composition_matrix,1);
num_lineages = size(subsampled_composition_matrix,2);
num_subjects = numel(unique({subsampled_composition_table.Subject{:}}));

mean_antagonism_frequency = mean(subsampled_ZOI_matrix,1);
mean_sensitivity_frequency = mean(subsampled_ZOI_matrix,2)';
mean_abundance = sum(subsampled_composition_matrix)./sum(subsampled_composition_matrix>0);

% Parse growth rate for each lineage in experiment
median_growth = [];
for n=1:num_lineages
    hold on
    lineage = lineage_labels(n);
    idxs = ([interactions_by_isolate.metadata.Lineage] == lineage);
    isolate_list = unique([interactions_by_isolate.metadata(idxs).Stock]);
    replicate_medians = [];
    replicate_stds = [];
    for i=1:numel(isolate_list)
        isolate = isolate_list(i);
        [replicate_medians(i),replicate_stds(i),growthrates,days] = get_lineage_growthrate(growth_table,isolate);
    end
    median_growth(n) = nanmedian(replicate_medians);
end

for n=1:num_lineages
    subject_idxs = find(subsampled_composition_matrix(:,n)>0);
    num_subjects_colonized(n) = numel(unique({subsampled_composition_table.Subject{subject_idxs}}));
end
prevalence = num_subjects_colonized./num_subjects;


% Expand vectors to treat samples independently for abundance predictors
abundance_vector = [];
antagonism_vector = [];
sensitivity_vector = [];
growth_vector = [];

for s=1:num_samples
    lineages_in_sample = find(subsampled_composition_matrix(s,:)>0); %nonzero
    abundance_vector = [abundance_vector, subsampled_composition_matrix(s,lineages_in_sample)];
    antagonism_vector = [antagonism_vector, mean_antagonism_frequency(lineages_in_sample)];
    sensitivity_vector = [sensitivity_vector, mean_sensitivity_frequency(lineages_in_sample)];
    growth_vector = [growth_vector,median_growth(lineages_in_sample)];
end
abundance_vector = abundance_vector';
antagonism_vector = antagonism_vector';
sensitivity_vector = sensitivity_vector';
growth_vector = growth_vector';


% Per lineage
correlation_structure.lineage_abundance = mean_abundance';
correlation_structure.lineage_antagonism = mean_antagonism_frequency';
correlation_structure.lineage_sensitivity = mean_sensitivity_frequency';
correlation_structure.lineage_growth = median_growth';
correlation_structure.lineage_prevalence = prevalence';
% Per Sample
correlation_structure.sample_abundance = abundance_vector;
correlation_structure.sample_antagonism = antagonism_vector;
correlation_structure.sample_sensitivity = sensitivity_vector;
correlation_structure.sample_growth = growth_vector;


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