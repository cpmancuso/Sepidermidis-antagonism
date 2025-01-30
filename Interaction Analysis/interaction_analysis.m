%% Interaction analysis script
% Chris Mancuso January 2025 
% This script is designed to produce all figures necessary for the analysis
% of interactions in "Intraspecies antagonism is a barrier to transmission
% in the skin microbiome". The script takes metadata, Zone of Inhibition
% measurements, and images as inputs. Proceed to the next code block to
% define parameters for this run. 

clc
clear all
close all
analysis_dir = pwd();
parent_dir = 'C:\Users\cmanc\Dropbox\Lieberman Lab\Personal lab notebooks\Chris Mancuso\Antagonism Paper\Revision Code';
addpath([analysis_dir '/myscripts'])
addpath([analysis_dir '/data_for_code'])

% Set default figure paramters
set(groot,'defaultAxesFontSize',8)
set(groot,'defaultAxesFontName','Helvetica')
set(groot,'defaultfigureposition',[100 100 350 350]) %approx 3.5in x 3.5in
set(groot,'defaultTextInterpreter','none')

% Point to data and load
disp('Loading data...')
folder_prefix = 'C:/Users/cmanc/Dropbox/Lieberman Lab/Personal lab notebooks/Chris Mancuso/Antagonism Paper/Publication Code/Image Analysis/';
folder_name={[folder_prefix 'AutoZOIs_median_filter - Block1'];[folder_prefix 'AutoZOIs_median_filter - Block1v2'];[folder_prefix 'AutoZOIs_median_filter - Block2v1'];[folder_prefix 'AutoZOIs_median_filter - Block2']};
datafiles={[folder_name{1} '/processed_data_2022_06_09_measured_16-Nov-2022.mat'];[folder_name{2} '/processed_data_2022_07_08_measured_16-Nov-2022.mat'];[folder_name{3} '/processed_data_2022_07_22_measured_16-Nov-2022.mat'];[folder_name{4} '/processed_data_2022_07_15_measured_16-Nov-2022.mat']};
image_dir={[folder_name{1} '/Masked Images/220606_T72_'];[folder_name{2} '/Masked Images/220705_T72_'];[folder_name{3} '/Masked Images/220719_T72_'];[folder_name{4} '/Masked Images/220712_T72_']};
interaction_structure = build_interaction_structure(datafiles);
metadata = readtable('updated_header_data.csv','TextType','string');

metadata.Name = strrep(metadata.Name,'_','-');

interaction_structure.metadata = table2struct(metadata); clear metadata
fnames = join(fieldnames(interaction_structure),', ');
disp(['Loaded ' fnames{:} ' into interaction structure.'])

%% Define parameters 
ZOI_depth_threshold = 8; %intensity depth, recommend 8
ZOI_noise_threshold = 4; %factor by which depth of ZOI must exceed noise (from isolated growth), recommend 4
ZOI_AUC_lower_threshold = 50; %minimum AUC value to be considered a true ZOI, recommend 50
ZOI_AUC_upper_threshold = 1000; %minimum AUC value to override depth & noise filters, recommend 1000
num_sims = 1000;

% set defaults, certain analyses may override
groupby_option = 'lineage'; %'isolate', 'lineage', 'none';
minmax_option = 'max'; %'min', 'max', or 'mode' interactions for group
remove_idiosyncratic = true; %choose whether to remove isolates with idiosyncratic sensitivity, recommend true
remove_superantagonists = false; %choose whether to remove isolates which antagonize most S. epi, recommend false
remove_replacements = false; % choose whether to remove isolates that were replaced for the M9 screen (absent from TSA screen), recommend false
combine_interaction_structure = false; % choose whether to combine interaction structures from multiple experiments, recommend false
combination_method = 'or'; % 'only1, only2, 'and', 'or', to choose how to combine interaction structures, not used if combine_interaction_structure is false

if ~remove_replacements
    interaction_structure = expand_interaction_structure(interaction_structure,'C:\Users\cmanc\Dropbox\Lieberman Lab\Personal lab notebooks\Chris Mancuso\Antagonism Paper\Revision Code\Image Analysis\replacement_TSA_plates_2024\automated_ZOI_results_TSA_extras_analyzed_26-Nov-2024.mat');
    for n=1:numel(interaction_structure.metadata)
        if strcmp(interaction_structure.metadata(n).Trustworthy,'replacement')
            interaction_structure.metadata(n).Trustworthy = "TRUE";
        end
    end
end

if remove_idiosyncratic
    for n=1:numel(interaction_structure.metadata)
        if ismember(interaction_structure.metadata(n).Stock,[211,313,292,332])
            interaction_structure.metadata(n).Trustworthy = 'idiosyncratic';
        end
    end
end

if remove_superantagonists
    for n=1:numel(interaction_structure.metadata)
        if ismember(interaction_structure.metadata(n).Lineage,[34,36,51,14])
            interaction_structure.metadata(n).Trustworthy = 'superantagonist';
        end
    end
end

%% Make ZOI calls then save interaction structure
% Make ZOI calls
disp('Making ZOI calls...')
interaction_structure.ZOI_call = make_ZOI_calls(interaction_structure,ZOI_depth_threshold,ZOI_noise_threshold,ZOI_AUC_lower_threshold,ZOI_AUC_upper_threshold);
% NOTE! 'interaction_structure' should remain unfiltered and sorted by index, create new variables for filtering, sorting purposes

save('interaction_structure_TSA.mat','interaction_structure')

%% Merge interaction structures from multiple analyses
if combine_interaction_structure
    additional_interactions = load('interaction_structure_M9.mat');
    load('lineage_order_expanded.mat')
    interaction_structure = merge_interaction_structure(interaction_structure,additional_interactions.interaction_structure,combination_method);
    if ~remove_replacements
        for n=1:numel(interaction_structure.metadata)
            if strcmp(interaction_structure.metadata(n).Trustworthy,'was replaced')
                interaction_structure.metadata(n).Trustworthy = "TRUE";
            end
        end
    end
else
    load('lineage_order_original.mat')
end

%% Group and sort interaction structures
% Group interactions and remove problematic samples
disp('Grouping interactions...')
% remove samples that were deemed untrustworthy for any reason (non "TRUE" values)
idxs = contains(vertcat(interaction_structure.metadata.Trustworthy),'TRUE');

subsample_structure = subsample_interaction_structure(interaction_structure,idxs);

% Create groupings by isolate, lineage, replicate
interactions_by_replicate = subsample_structure;
interactions_by_isolate = group_interactions(subsample_structure,'isolate',minmax_option);
interactions_by_lineage = group_interactions(interactions_by_isolate,'lineage',minmax_option);
clear subsample_structure

% Sort Interaction Structure
interactions_by_replicate = sort_interaction_structure(interactions_by_replicate,{'Phylogroup'});
interactions_by_isolate = sort_interaction_structure(interactions_by_isolate,{'Phylogroup'});
interactions_by_lineage = sort_interaction_structure(interactions_by_lineage,{'Lineage'},lineage_order);

idxs = contains(vertcat(interactions_by_lineage.metadata.Species),'epidermidis');
interactions_by_lineage_sepi = subsample_interaction_structure(interactions_by_lineage,idxs);
interactions_by_lineage_sepi = sort_interaction_structure(interactions_by_lineage_sepi,{'Lineage'},lineage_order);

idxs = contains(vertcat(interactions_by_isolate.metadata.Species),'epidermidis');
interactions_by_isolate_sepi = subsample_interaction_structure(interactions_by_isolate,idxs);
interactions_by_isolate_sepi = sort_interaction_structure(interactions_by_isolate_sepi,{'Lineage'},lineage_order);

interactions_by_isolate = sort_interaction_structure(interactions_by_isolate,{'Lineage','Name'},lineage_order);

%% FIGURE GENERATION %%
% Each section below is organized by figure panel. Set fignums to distinct
% values for each figure you'd like to generate. The following information are
% provided for each section

% Generates: the figure this data corresponds to in the paper
% Description: high level description of the analysis being performed
% Requirements: explanation for any specific varaiable requirements that would causes errors
%

%% Repopulate composition table and plot relative abundance
% Generates: Supplementary Figure 2
% Description: Plots relative lineage abundance for all subjects in the study
% Requirements: None
%

fignum=0; % Makes one figure per family
if fignum
    read_and_plot_composition('all_isolate_lineages.csv','updated_header_data.csv','Sepi_lineage_level_frequencies_filtered.csv',false,fignum)
end
% load composition table once this is done, or if skipped
load('composition_table.mat')


%% Make clickable interaction heatmap
% Generates: None, see heatmap with overlay below
% Description: Makes clickable heatmap for data exploration where you can bring up specific spot images
% Requirements: set click_option to choose between lineage and isolate levels
%
fignum = 0;
if fignum
    click_option = 'lineage'; %lineage or isolate
    switch click_option
        case 'isolate'
            fighandle = clickable_interaction_heatmap(interactions_by_isolate,interaction_structure,'Name',image_dir,fignum,ZOI_depth_threshold);
        case 'lineage'
            fighandle = clickable_interaction_heatmap(interactions_by_lineage,interaction_structure,'Lineage',image_dir,fignum,ZOI_depth_threshold);
    end
end

%% Make heatmap with overlay
% Generates: Figure 2A, Supplementary Figure 3,8
% Description: Makes interaction heatmap with red overlays at specified level
% Requirements: set overlay_option to choose between lineage and isolate levels
%

fignum = 0;

if fignum
    overlay_option = 'isolate'; %lineage or isolate
    switch overlay_option
        case 'lineage'
            fighandle = overlay_interaction_heatmap(interactions_by_lineage,'Lineage',fignum,'Phylogroup');
        case 'isolate'
            fighandle = overlay_interaction_heatmap(interactions_by_isolate,'Name',fignum,'Species');
            export_interaction_table(interactions_by_isolate,'Name','interactions_by_isolate.csv');
    end
end

%% Compare heatmaps between two experiments
% Generates: Supplementary Figure 9
% Description: Plot heatmap showing interaction differences between two experiments
% Requirements: requires combined interaction structure to be defined, see above
%

fignum = 0;
if fignum
    if ~combine_interaction_structure
       error('Error: Cannot compare experiments without combined interaction structure, set combine_interaction_structure...')
    end
    fighandle = compare_interaction_heatmaps(interactions_by_lineage_sepi,'Lineage',fignum,'Species');
end

%% Plot composition by antagonism
% Generates: Supplementary Figure 12
% Description: Plot composition of each family, colored by amount of antagonism
% Requirements: Cannot currently handle partial data from replacements
%

fignum=0; %makes one figure per family, 
if fignum 
    if ~remove_replacements
       error('Error: partial interaction data not supported for this analysis, set remove_replacements...')
    end
    plot_composition_by_antagonism(composition_table,interactions_by_lineage_sepi,fignum)
end


%% Make ANI plots
% Generates: Supplementary Figure 5
% Description: Plot antagonism vs ANI and AAI
% Requirements: Cannot currently handle partial data from replacements
%

fignum = 0;
if fignum
    % subsample only isolates with fastANI data
    idxs = (vertcat(interactions_by_lineage.metadata.AAI_index)>0)|(vertcat(interactions_by_lineage.metadata.ANI_index)>0);
    interactions_with_identity = subsample_interaction_structure(interactions_by_lineage,idxs);
    [interaction_structure_with_identity_out,fighandle] = plot_AAI_ANI(interactions_with_identity,'AAI_matrix.mat','fastANI_matrix.mat', fignum);
end

%% Compare targets of antagonism
% Generates: Figure 2B,C,D, Figure 3B, Supplementary Figure 5G
% Description: Plot antagonism vs ANI and AAI
% Requirements: Cannot currently handle partial data from replacements
%
% Figure 2B,C,D
shuffle_option = 'both'; % shuffle_option: rows, cols, or both to shuffle antagonizers and/or baits

% Phylogroup
fignum = 0;
if fignum
    % [fighandle p fishertable] = plot_fishers_exact(interactions_by_lineage_sepi,'Phylogroup','Phylogroup',fignum);
    [simulation_structure] = permute_interactions_by_group(interactions_by_lineage_sepi,'Phylogroup',shuffle_option,1000,fignum);
    [fighandle] = plot_between_group_differences(simulation_structure,'Phylogroup','Phylogroup',fignum);
end

% Agr Type
fignum = 0;
if fignum
    % [fighandle p fishertable] = plot_fishers_exact(interactions_by_lineage_sepi,'Agr_Type','agr Type',fignum);
    [simulation_structure] = permute_interactions_by_group(interactions_by_lineage_sepi,'Agr_Type',shuffle_option,1000,fignum);
    [fighandle] = plot_between_group_differences(simulation_structure,'Agr_Type','agr Type',fignum);
end
%
% Species
fignum = 0;
if fignum
    % [fighandle p fishertable] = plot_fishers_exact(interactions_by_lineage,'Species','Species',fignum);
    [simulation_structure] = permute_interactions_by_group(interactions_by_lineage,'Species',shuffle_option,1000,fignum);
    [fighandle] = plot_between_group_differences(simulation_structure,'Species','Species',fignum);
end

% Sample / Timepoint
fignum = 0;
if fignum
    [simulation_structure] = permute_interactions_by_timepoint(interactions_by_lineage_sepi,composition_table,shuffle_option,1000);
    [fighandle] = plot_between_group_differences(simulation_structure,'Sample','Sample',fignum);
end

% Lineage
fignum = 0;
if fignum
    [simulation_structure] = permute_interactions_by_group(interactions_by_isolate_sepi,'Lineage',shuffle_option,1000,fignum);
    [fighandle] = plot_between_group_differences(simulation_structure,'Lineage','Lineage',fignum);
end

%% Calculate antagonism frequency for each subject
% Generates: Figure 3A, Supplementary Figure 10A,B,C
% Description: Plot per sample Antagonism Frequency (AF) or mean AF per subject, family etc
% Requirements: None, supports partial data from replacements
% Note: p-values reported from t-tests, but assumptions are not valid

% remove composition data for unrepresented lineages
interactions_by_lineage = sort_interaction_structure(interactions_by_lineage,{'Lineage'},lineage_order);
idxs = contains(vertcat(interactions_by_lineage.metadata.Species),'epidermidis');
interactions_by_lineage_sepi = subsample_interaction_structure(interactions_by_lineage,idxs);
rep_lineages = [interactions_by_lineage_sepi.metadata.Lineage];
composition_matrix = composition_table{:,rep_lineages+3}; %filter out unrepresented lineages
ZOI_matrix = interactions_by_lineage_sepi.ZOI_call;

[weighted_freq_structure] = calculate_interaction_frequency(composition_matrix,ZOI_matrix,composition_table.Subject,'weighted');
[nonweighted_freq_structure] = calculate_interaction_frequency(composition_matrix,ZOI_matrix,composition_table.Subject,'nonweighted');

fignum = 0;
if fignum
    % Plot per sample
    [fighandle,p] = plot_interaction_frequency_stem(weighted_freq_structure.expected_freq,weighted_freq_structure.per_sample_interaction_freq,weighted_freq_structure.per_sample_interaction_freq_upperbound,composition_table.SID,'Samples',fignum);
    % [fighandle,pnw] = plot_interaction_frequency_stem(nonweighted_freq_structure.expected_freq,nonweighted_freq_structure.per_sample_interaction_freq,nonweighted_freq_structure.per_sample_interaction_freq_upperbound,composition_table.SID,'Samples',fignum+10);
end

fignum = 0;
if fignum
    % Plot per subject
    [fighandle,p] = plot_interaction_frequency_stem(weighted_freq_structure.expected_freq,weighted_freq_structure.per_subject_interaction_freq,[],weighted_freq_structure.subjects,'Subjects',fignum);
    % [fighandle,pnw] = plot_interaction_frequency_stem(nonweighted_freq_structure.expected_freq,nonweighted_freq_structure.per_subject_interaction_freq,[],nonweighted_freq_structure.subjects,'Subjects',fignum+10);

end

fignum = 0;
if fignum
    families = {'1','2','4','5','7','8'};
    for f=1:numel(families)
        per_family_interaction_freq(f) = mean(weighted_freq_structure.per_subject_interaction_freq(contains(weighted_freq_structure.subjects,families{f})));
        per_family_interaction_freq_nw(f) = mean(nonweighted_freq_structure.per_subject_interaction_freq(contains(nonweighted_freq_structure.subjects,families{f})));
    end
    [fighandle,p] = plot_interaction_frequency_stem(weighted_freq_structure.expected_freq,per_family_interaction_freq,[],families,'Families',fignum);
    % [fighandle,pnw] = plot_interaction_frequency_stem(nonweighted_freq_structure.expected_freq,per_family_interaction_freq_nw,[],families,'Families',fignum+10);
end

%% Permutation tests, shuffle across all families
% Generates: Figure 3C, Supplementary Figure 9B-E, 10D-G, 13
% Description: plots deltaAF observed vs simulations for 4 permutation analyses with different assumptions, 
% Requirements: Cannot currently handle partial data from replacements
%
weight_option = 'weighted'; % 'weighted' or 'nonweighted'; consider relative abundance? recommend weighted
replace_option = 'noreplace'; % 'replace' or 'noreplace'; shuflling with or without replacement, recommend noreplace

% Optionally shuffle interaction matrix before analysis
randomize_interactions = false; % strongly recommend false
if randomize_interactions
    disp('Randomized interactions...')
    interactions_for_permutation = interactions_by_lineage_sepi;
    interactions_for_permutation.ZOI_call = interactions_for_permutation.ZOI_call(reshape(randperm(numel(interactions_for_permutation.ZOI_call)),size(interactions_for_permutation.ZOI_call)));
else
    interactions_for_permutation = interactions_by_lineage_sepi;
end

% Shuffle lineages across all families, by subject (breaking composition structure)
% Supplementary Figure 10E 
fignum = 0;
if fignum
    figure(fignum)
    clf(fignum)
    if ~remove_replacements
       error('Error: partial interaction data not supported for this analysis, set remove_replacements...')
    end
    % remove composition data for unrepresented lineages
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
    [weighted_freq_structure] = calculate_interaction_frequency(subsampled_composition_matrix,subsampled_ZOI_matrix,subsampled_composition_table.Subject,weight_option);

    % permute, may take a few seconds
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_for_permutation, comp_idxs);
    [weighted_sims_structure] = permute_composition(subsampled_composition_table,subsampled_ZOI_matrix,num_sims,weight_option,replace_option,'population','subject');
    [fighandle,deltaAF_structure, pvals] = plot_simulation_results(weighted_freq_structure,weighted_sims_structure,num_sims,fignum);
end

% Shuffle lineages across all families, by row (retaining composition structure)
% Supplementary Figure 10D, 9B-E
fignum = 0;
if fignum
    figure(fignum)
    clf(fignum)
    if ~remove_replacements
       error('Error: partial interaction data not supported for this analysis, set remove_replacements...')
    end
    % remove composition data for unrepresented lineages
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
    [weighted_freq_structure] = calculate_interaction_frequency(subsampled_composition_matrix,subsampled_ZOI_matrix,subsampled_composition_table.Subject,weight_option);

    % permute, may take a few seconds
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_for_permutation, comp_idxs);
    [weighted_sims_structure] = permute_composition(subsampled_composition_table,subsampled_ZOI_matrix,num_sims,weight_option,replace_option,'population','column');
    [fighandle,deltaAF_structure, pvals] = plot_simulation_results(weighted_freq_structure,weighted_sims_structure,num_sims,fignum);
    
end

% Shuffle lineages within families, by subject (breaking composition structure)
% Supplementary Figure 10G 
fignum = 0;
if fignum
    figure(fignum)
    clf(fignum)
    if ~remove_replacements
       error('Error: partial interaction data not supported for this analysis, set remove_replacements...')
    end
    % remove composition data for unrepresented lineages
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
    [weighted_freq_structure] = calculate_interaction_frequency(subsampled_composition_matrix,subsampled_ZOI_matrix,subsampled_composition_table.Subject,weight_option);

    % permute, may take a few seconds
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_for_permutation, comp_idxs);
    [weighted_sims_structure] = permute_composition(subsampled_composition_table,subsampled_ZOI_matrix,num_sims,weight_option,replace_option,'family','subject');
    [fighandle,deltaAF_structure, pvals] = plot_simulation_results(weighted_freq_structure,weighted_sims_structure,num_sims,fignum);
end

% Shuffle lineages within families, by row (retaining composition structure)
% Supplementary Figure 10F 
fignum = 0;
if fignum
    figure(fignum)
    clf(fignum)
    if ~remove_replacements
       error('Error: partial interaction data not supported for this analysis, set remove_replacements...')
    end
    % remove composition data for unrepresented lineages
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
    [weighted_freq_structure] = calculate_interaction_frequency(subsampled_composition_matrix,subsampled_ZOI_matrix,subsampled_composition_table.Subject,weight_option);

    % permute, may take a few seconds
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_for_permutation, comp_idxs);
    [weighted_sims_structure] = permute_composition(subsampled_composition_table,subsampled_ZOI_matrix,num_sims,weight_option,replace_option,'family','column');
    [fighandle,deltaAF_structure, pvals] = plot_simulation_results(weighted_freq_structure,weighted_sims_structure,num_sims,fignum);    
end


%% Make subject-ordinated heatmaps to compare families
% Generates: Figure 3D, Supplementary Figure 11
% Description: plots heatmaps with duplicated lineages / isolates for each sample included in list
% Requirements: set subject_option to choose which subjects to include
%

fignum=0;
if fignum
    subject_option = '57';
    % create duplicated ZOI data for the following subjects
    switch subject_option
        case '57'
            samples_to_include = {'5PA1','5PB3','7AA1','7AB4','7PA1'};
        case 'all'
            samples_to_include = {'1AA1','1AA3','1AA4','1PA1','1PA2','1PA4','1PB2','1PB3','2AA1','2AA3','2PA1','2PA3','2PB1','4AA1','4AB1','5PA1','5PB3','7AA1','7AA4','7AB4','7PA1','8AA4','8AB4','8AC4','8PA3','8PB1','8PB4'};
    end
    included_idxs = [];
    included_labels = {};
    replicated_subject_labels = {};
    num = 0;
    for s = 1:numel(samples_to_include)
        temp_composition_matrix = composition_table{(find(contains(composition_table.SID,samples_to_include{s}))),4:end};
        lineages_represented = find(temp_composition_matrix(1,:)>0);
        for l=1:numel(lineages_represented)
            lineage = lineages_represented(l);
            idx = find([interactions_by_lineage_sepi.metadata.Lineage]==lineage);
            if ~isempty(idx)
                num = num+1;
                included_idxs(num) = idx;
                included_labels{num} = [samples_to_include{s} ' - ' num2str(lineage)];
                replicated_subject_labels{num} = samples_to_include{s};
            end
        end
    end
    % make new interaction structure with replicated lineages
    interactions_with_replication = subsample_interaction_structure(interactions_by_lineage,included_idxs);
    for n=1:num
        interactions_with_replication.metadata(n).Name = included_labels{n};
        interactions_with_replication.metadata(n).Subject = replicated_subject_labels{n};
    end
    fighandle = overlay_interaction_heatmap(interactions_with_replication,'Name',fignum,'Subject');
end

%% Make heatmap to show idiosyncratic isolates
% Generates: Figure 4A
% Description: plots heatmap with only specified lineages, in this case, idiosyncratic lineages 20, 37, and 58
% Requirements: set lineages_to_include
%

fignum=0;
if fignum
    lineages_to_include = [7,12,21,35,57,20,37,58];
    idiosyncratic_idxs = [];
    num = 0;
    for l = 1:numel(lineages_to_include)
        idxs = find([interactions_by_isolate_sepi.metadata.Lineage] == lineages_to_include(l));
        idiosyncratic_idxs = [idiosyncratic_idxs idxs];
    end
    interactions_idiosyncratic = subsample_interaction_structure(interactions_by_isolate_sepi,idiosyncratic_idxs);
    fighandle = overlay_interaction_heatmap(interactions_idiosyncratic,'Name',fignum,'Lineage');
end

%% Linear and spearman correlations
% Generates: Figure 3E,F, Supplementary Figure 15 C,D,E
% Description: plots linear and spearman correlations for sets of data, see each
% Requirements: Cannot currently handle partial data from replacements
%

% Generate structure with data to be correlated
% remove composition data for unrepresented lineages
growth_filename = "figure_S12_growth_rates_2308_2312";
comp_idxs = 1:size(composition_table,1);
[subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix,lineage_labels] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
correlation_structure = prepare_for_correlations(subsampled_composition_table,subsampled_ZOI_matrix,growth_filename,interactions_by_isolate,lineage_labels);

% Sharing vs antagonism
% Figure 3E
fignum = 0; 
if fignum
    [fighandle,max_ants_per_sample] = plot_sharing_vs_antagonism(subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix,fignum);
end

% Abundance vs antagonism
% Figure 3F
fignum = 0; 
if fignum
    [fighandle] = plot_abundance_vs_antagonism(correlation_structure,lineage_labels,fignum);
end

% Abundance vs growth rate
% Supplementary Figure 15C
fignum = 0; 
if fignum
    [fighandle] = plot_abundance_vs_growthrate(correlation_structure,lineage_labels,fignum);
end

% Prevalence vs antagonism
% Supplementary Figure 15D
fignum = 0; 
if fignum
    [fighandle] = plot_prevalence_vs_antagonism(correlation_structure,lineage_labels,fignum);
end

% Growth rate vs antagonism
% Supplementary Figure 15E
fignum = 0; 
if fignum
    fighandle = plot_growthrate_vs_antagonism(correlation_structure,lineage_labels,fignum);
end

%% Plot fold enrichment of different gene categories
% Generates: Supplementary Figure 1
% Description: plots fold enrichment of COGs / gene lists
% Requirements: None

fignum=0; %makes 2 figs
if fignum
    accessory_stats = readtable("figure_S1_accessory_stats_output.csv");
    gainloss_stats = readtable("figure_S1_gainloss_stats_output.csv");
    fighandle = plot_core_fold_enrichment(fignum,accessory_stats);
    fighandle = plot_core_fold_enrichment(fignum+1,gainloss_stats);
end

%% Cluster interactions by similarity and plot dendrograms
% Generates: Supplementary Figure 6
% Description: plots clustergram of interactions and lineages
% Requirements: None, but version in manuscript included idiosyncratic isolates

fignum=0; % clustergram is buggy, opens multiple figs
if fignum
    if remove_idiosyncratic
       warning('Warning: Manuscript interactions were clustered with remove_idiosyncratic = false...')
    end
    fighandle = plot_interaction_clustergram(interactions_by_lineage,fignum);
end

%% Mechanism upset plot
% Generates: Supplementary Figure 7
% Description: creates upset plot of preliminary mechanism screen data
% Requirements: None

fignum=0;
if fignum
    fighandle = plot_mechanism_upset('mechanism_upset_data.csv',fignum);
end

%% Plot composition by agr type and phylogroup
% Generates: Supplementary Figure 14
% Description: plots composition of each sample by phylogroup and by agr type
% Requirements: None

fignum=0;
if fignum
    fighandle = plot_agr_and_phylo_composition(interactions_by_lineage_sepi,rep_lineages,composition_table,composition_matrix,fignum);
end

%% Plot MIC for relevant isolates
% Generates: Figure 4C, Supplementary Figure 17
% Description: plots MIC for selected compounds against selected isolates  
% Requirements: None

% Figure 4
fignum = 0;
if fignum
    fighandle = plot_MIC_change('figure_4_MICs_with_revertant.xlsx',fignum);
end

% Supplementary Figure 17
fignum=0;
if fignum
    fighandle = summarize_MIC('figure_S11_MICs_summary.xlsx',fignum);
end

%% Plot growth rate by lineage and of selected isolates
% Generates: Supplementary Figure 15D,E
% Description: plots growth rate of selected isolates  
% Requirements: None

% Supplementary Figure 15E
fignum = 0;
if fignum
    growth_filename = "figure_S12_idiosyncratic_growth";
    fighandle = plot_idiosyncratic_growth(growth_filename,fignum);
end


%% Look for cases of intralineage variation
% Generates: Supplementary Figure 16A
% Description: plots intralineage differences in antagonism or sensitivity and generates a table with significant hits  
% Requirements: remove_idiosyncratic must be false in order to include all isolates, cannot handle partial data from replacements
fignum=0;
if fignum
    if remove_idiosyncratic
       error('Error: Must include all isolates, set remove_idiosyncratic = false...')
    end
    if ~remove_replacements
       error('Error: partial interaction data not supported for this analysis, set remove_replacements...')
    end
    [fighandle, ant_outliers_table, sen_outliers_table] = plot_intralineage_variation(interactions_by_isolate,fignum);
end

%% Summarize findings
% Generates: None
% Description: Displays useful statistics  
% Requirements: none

summarize_option = true;
if summarize_option
    % Ratio of shared / unshared lineages
    sharing_table = grpstats(composition_matrix,composition_table{:,3},"mean");
    num_shared = sum(sum(sharing_table>0)>1);
    num_unshared = sum(sum(sharing_table>0)>0)-num_shared;

    % Superantagonizers vs nonanatagonizers
    num_ants = nansum(ZOI_matrix);
    is_superant = num_ants>(numel(num_ants)./2);
    superant_fraction = nansum(num_ants(is_superant))./sum(num_ants);
    nonant_fraction = nansum(num_ants==0)./numel(num_ants);
    disp([num2str(sum(is_superant),3) ' superantagonist lineages are responsible for ' num2str(superant_fraction,3) ' of antagonisms.'])
    disp([num2str(nonant_fraction,3) ' of lineages do not antagonize any others.']) 
    disp([num2str(sum(num_ants)./(numel(ZOI_matrix)-size(ZOI_matrix,1)),3) ' of interactions are antagonistic.'])
    disp(['Weighted Antagonism Frequency across the cohort is ' num2str(weighted_freq_structure.expected_freq,3)])  
end

%% Fit antagonism data to distribution
% Generates: Supplementary Figure 13
% Description: plots histogram of antagonism by lineages and fits exponential distribution  
% Requirements: None

fignum=0;
if fignum
    [fighandle,hExp, pExp] = plot_exponential_dist(num_ants,fignum);
end

%% Other useful analyses
% Other analyses which may be useful but were not included in the paper

% Do children have higher AF? No
fignum=0;
if fignum
    freq = weighted_freq_structure.per_subject_interaction_freq;
    AF_child = freq(~contains(weighted_freq_structure.subjects,'P'));
    AF_adult = freq(contains(weighted_freq_structure.subjects,'P'));
    [H,P] = ttest2(AF_child,AF_adult)
end

% Do antagonizers affect shannon diversity? No
fignum = 0;
if fignum
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
    fighandle = plot_shannon_vs_antagonism(subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix,fignum);
end

% Compare dMRCA for antagonists vs non-antagonists. Not significant
fignum = 0;
if fignum
    fighandle = plot_dMRCA_vs_antagonism(interactions_by_lineage_sepi,composition_table,'dMRCA_table.csv','antagonism',fignum);
    fighandle = plot_dMRCA_vs_antagonism(interactions_by_lineage_sepi,composition_table,'dMRCA_table.csv','sensitivity',fignum+1);
end
