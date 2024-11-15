%% Interaction analysis script
% Chris Mancuso April 2024 
% This script is designed to produce all figures necessary for the analysis
% of interactions in "Intraspecies antagonism is a barrier to transmission
% in the skin microbiome". The script takes metadata, Zone of Inhibition
% measurements, and images as inputs. Proceed to the next code block to
% define parameters for this run. Set fignum to different values in order
% to generate figures.

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
folder_prefix = [parent_dir, '/Image Analysis/'];
folder_name={[folder_prefix 'Block_1v1_241104'];[folder_prefix 'Block_1v2_241106'];[folder_prefix 'Block_2v1_241106'];[folder_prefix 'Block_2v2_241104']};
datafiles={[folder_name{1} '/automated_ZOI_results_2024_10_20_block1v1_analyzed_07-Nov-2024.mat'];[folder_name{2} '/automated_ZOI_results_2024_10_26_block1v2_analyzed_08-Nov-2024.mat'];[folder_name{3} '/automated_ZOI_results_2024_11_02_block2v1_analyzed_08-Nov-2024.mat'];[folder_name{4} '/automated_ZOI_results_2024_10_18_block2v2_analyzed_08-Nov-2024.mat']};
image_dir={[folder_name{1} '/Masked Images/Block1v1_24_10_20_'];[folder_name{2} '/Masked Images/Block1v2_24_10_26_'];[folder_name{3} '/Masked Images/Block2v1_24_11_02_'];[folder_name{4} '/Masked Images/Block2v2_24_10_18_']};
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
minmax_option = 'max'; %'min' or 'max' interactions for group
collapse_non_epi = false; %choose whether to collapse non S. epidermidis isolates into a single "lineage" per species, recommend false
remove_idiosyncratic = false; %choose whether to remove isolates with idiosyncratic sensitivity, recommend false
remove_superantagonists = false; %choose whether to remove isolates which antagonize most S. epi, recommend false
remove_replacements = false; % choose whether to remove isolates that were replaced for the M9 screen (absent from TSA screen)
if collapse_non_epi
    for n=1:numel(interaction_structure.metadata)
        interaction_structure.metadata(n).Lineage = round(interaction_structure.metadata(n).Lineage);
    end
end

if remove_idiosyncratic
    for n=1:numel(interaction_structure.metadata)
        if ismember(interaction_structure.metadata(n).Stock,[211,313,227,311])
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

if ~remove_replacements
    for n=1:numel(interaction_structure.metadata)
        if strcmp(interaction_structure.metadata(n).Trustworthy,'replacement')
            interaction_structure.metadata(n).Trustworthy = 'TRUE';
        end
    end
end

%% Make ZOI calls then save interaction structure
% Make ZOI calls
disp('Making ZOI calls...')
interaction_structure.ZOI_call = make_ZOI_calls(interaction_structure,ZOI_depth_threshold,ZOI_noise_threshold,ZOI_AUC_lower_threshold,ZOI_AUC_upper_threshold);
% NOTE! 'interaction_structure' should remain unfiltered and sorted by index, create new variables for filtering, sorting purposes

save('interaction_structure_M9.mat','interaction_structure')

%% Optional: Merge interaction structures from multiple analyses
merge_interactions = False;
if merge_interactions
    additional_interactions = load('interaction_structure_TSA.mat');
    interaction_structure = merge_interaction_structure(interaction_structure,additional_interactions.interaction_structure,'or');
end

%% Group interactions and remove problematic samples
disp('Grouping interactions...')
% remove samples that were deemed untrustworthy for any reason (non "TRUE" values)
idxs = contains(vertcat(interaction_structure.metadata.Trustworthy),'TRUE');

subsample_structure = subsample_interaction_structure(interaction_structure,idxs);

% Create groupings by isolate, lineage, replicate
interactions_by_replicate = subsample_structure;
interactions_by_isolate = group_interactions(subsample_structure,'isolate',minmax_option);
interactions_by_lineage = group_interactions(subsample_structure,'lineage',minmax_option);
clear subsample_structure

%% Sort Interaction Structure
interactions_by_replicate = sort_interaction_structure(interactions_by_replicate,{'Phylogroup'});
interactions_by_isolate = sort_interaction_structure(interactions_by_isolate,{'Phylogroup'});
interactions_by_lineage = sort_interaction_structure(interactions_by_lineage,{'Phylogroup'});


%% Repopulate composition table and plot relative abundance
% Supplementary Figure 2
fignum=0; %makes several figures
if fignum
    read_and_plot_composition('all_isolate_lineages.csv','all_header_data.csv','Sepi_lineage_level_frequencies_filtered.csv',false,fignum)
end
% load composition table once this is done, or if skipped
load('composition_table.mat')


%% Make clickable interaction heatmap
% For data exploration

fignum = 0;
if fignum
    fighandle = clickable_interaction_heatmap(interactions_by_isolate,interaction_structure,'Name',image_dir,fignum,ZOI_depth_threshold);
    % fighandle = clickable_interaction_heatmap(interactions_by_lineage,interaction_structure,'Lineage',image_dir,fignum,ZOI_depth_threshold);
end
%% Make heatmap with overlay
% To generate Figure 2A, Supplementary Figure 3 
fignum = 0;

idxs = contains(vertcat(interactions_by_lineage.metadata.Species),'epidermidis');
interactions_by_lineage_sepi = subsample_interaction_structure(interactions_by_lineage,idxs);
interactions_by_lineage_sepi = sort_interaction_structure(interactions_by_lineage_sepi,{'Lineage'});

idxs = contains(vertcat(interactions_by_isolate.metadata.Species),'epidermidis');
interactions_by_isolate_sepi = subsample_interaction_structure(interactions_by_isolate,idxs);
interactions_by_isolate_sepi = sort_interaction_structure(interactions_by_isolate_sepi,{'Lineage'});

interactions_by_isolate = sort_interaction_structure(interactions_by_isolate,{'Lineage','Name'});


if fignum
    fighandle = overlay_interaction_heatmap(interactions_by_lineage,'Name',fignum,'Species');
end

%% Plot composition by antagonism
% Unused analysis
fignum=0; %makes several figures
if fignum
    plot_composition_vs_antagonism(composition_table,interactions_by_lineage_sepi,fignum)
end

%% Make ANI plots
% To generate Supplementary Figure 5

% subsample only isolates with fastANI data
idxs = (vertcat(interactions_by_lineage.metadata.AAI_index)>0)|(vertcat(interactions_by_lineage.metadata.ANI_index)>0);
interactions_with_identity = subsample_interaction_structure(interactions_by_lineage,idxs);

fignum = 0;
if fignum
    [interaction_structure_with_identity_out,fighandle] = plot_AAI_ANI(interactions_with_identity,'AAI_matrix.mat','fastANI_matrix.mat', fignum);
end
%% Calculate interaction frequency for each subject
% To generate Figure 3A, and Supplementary Figure 9A,B,C

% remove composition data for unrepresented lineages
interactions_by_lineage = sort_interaction_structure(interactions_by_lineage,{'Lineage'});
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
    figure(fignum)
    clf(fignum)
    [fighandle,p] = plot_interaction_frequency_stem(weighted_freq_structure.expected_freq,weighted_freq_structure.per_sample_interaction_freq,composition_table.SID,'Samples',fignum);
    [fighandle,pnw] = plot_interaction_frequency_stem(nonweighted_freq_structure.expected_freq,nonweighted_freq_structure.per_sample_interaction_freq,composition_table.SID,'Samples',fignum+10);
    disp(['Per sample Weighted p=' num2str(p) ' Non-Weighted p=' num2str(pnw)])
end

fignum = 0;
if fignum
    % Plot per subject
    [fighandle,p] = plot_interaction_frequency_stem(weighted_freq_structure.expected_freq,weighted_freq_structure.per_subject_interaction_freq,weighted_freq_structure.subjects,'Subjects',fignum);
    [fighandle,pnw] = plot_interaction_frequency_stem(nonweighted_freq_structure.expected_freq,nonweighted_freq_structure.per_subject_interaction_freq,nonweighted_freq_structure.subjects,'Subjects',fignum+10);
    disp(['Per subject Weighted p=' num2str(p) ' Non-Weighted p=' num2str(pnw)])

end

fignum = 0;
if fignum
    figure(fignum)
    clf(fignum)
    families = {'1','2','4','5','7','8'};
    for f=1:numel(families)
        per_family_interaction_freq(f) = mean(weighted_freq_structure.per_subject_interaction_freq(contains(weighted_freq_structure.subjects,families{f})));
        per_family_interaction_freq_nw(f) = mean(nonweighted_freq_structure.per_subject_interaction_freq(contains(nonweighted_freq_structure.subjects,families{f})));
    end
    figure(fignum)
    clf(fignum)
    [fighandle,p] = plot_interaction_frequency_stem(weighted_freq_structure.expected_freq,per_family_interaction_freq,families,'Families',fignum);
    [fighandle,pnw] = plot_interaction_frequency_stem(nonweighted_freq_structure.expected_freq,per_family_interaction_freq_nw,families,'Families',fignum+10);
    disp(['Per family Weighted p=' num2str(p) ' Non-Weighted p=' num2str(pnw)])
end



%% Calculate interaction frequency for each family separately
% Unused analysis

idxs = contains(vertcat(interactions_by_lineage.metadata.Species),'epidermidis');
interactions_by_lineage_sepi = subsample_interaction_structure(interactions_by_lineage,idxs);
fignum = 0;
if fignum
    fighandle = plot_interaction_frequency_for_each_family(composition_table,interactions_by_lineage,fignum);
end


%% Permutation tests, shuffle across all families and plot CDFs
% Four analyses with different assumptions, see figure for each

weight_option = 'weighted'; % 'weighted' or 'nonweighted'; consider relative abundance? recommend weighted
replace_option = 'noreplace'; % 'replace' or 'noreplace'; shuflling with or without replacement, recommend noreplace

% Optionally shuffle interaction matrix before analysis
randomize_interactions = false; % recommend false
if randomize_interactions
    interactions_for_permutation = interactions_by_lineage_sepi;
    interactions_for_permutation.ZOI_call = interactions_for_permutation.ZOI_call(reshape(randperm(numel(interactions_for_permutation.ZOI_call)),size(interactions_for_permutation.ZOI_call)));
else
    interactions_for_permutation = interactions_by_lineage_sepi;
end

% Shuffle lineages across all families, by row (retaining composition structure)
% Supplementary Figure 9D 
fignum = 0;
if fignum
    figure(fignum)
    clf(fignum)
    % remove composition data for unrepresented lineages
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_for_permutation, comp_idxs);
    
    % permute, may take a few seconds
    [weighted_freq_structure] = calculate_interaction_frequency(subsampled_composition_matrix,subsampled_ZOI_matrix,subsampled_composition_table.Subject,weight_option);
    [weighted_sims_structure] = permute_composition(subsampled_composition_table,subsampled_ZOI_matrix,num_sims,weight_option,replace_option,'population','column');
    [fighandle,deltaIF_structure, pvals] = plot_simulation_results(weighted_freq_structure,weighted_sims_structure,num_sims,fignum);
    
end

% Shuffle lineages across all families, by element (breaking composition structure)
% Supplementary Figure 9E 
fignum = 0;
if fignum
    figure(fignum)
    clf(fignum)
    % remove composition data for unrepresented lineages
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_for_permutation, comp_idxs);
    
    % permute, may take a few seconds
    [weighted_freq_structure] = calculate_interaction_frequency(subsampled_composition_matrix,subsampled_ZOI_matrix,subsampled_composition_table.Subject,weight_option);
    [weighted_sims_structure] = permute_composition(subsampled_composition_table,subsampled_ZOI_matrix,num_sims,weight_option,replace_option,'population','subject');
    [fighandle,deltaIF_structure, pvals] = plot_simulation_results(weighted_freq_structure,weighted_sims_structure,num_sims,fignum);
end


% Shuffle lineages within families, by row (retaining composition structure)
% Supplementary Figure 9F 
fignum = 0;
if fignum
    figure(fignum)
    clf(fignum)
    % remove composition data for unrepresented lineages
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_for_permutation, comp_idxs);
    
    % permute, may take a few seconds
    [weighted_freq_structure] = calculate_interaction_frequency(subsampled_composition_matrix,subsampled_ZOI_matrix,subsampled_composition_table.Subject,weight_option);
    [weighted_sims_structure] = permute_composition(subsampled_composition_table,subsampled_ZOI_matrix,num_sims,weight_option,replace_option,'family','column');
    [fighandle,deltaIF_structure, pvals] = plot_simulation_results(weighted_freq_structure,weighted_sims_structure,num_sims,fignum);    
end


% Shuffle lineages within families, by element (breaking composition structure)
% Supplementary Figure 9G 
fignum = 0;
if fignum
    figure(fignum)
    clf(fignum)
    % remove composition data for unrepresented lineages
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_for_permutation, comp_idxs);
    
    % permute, may take a few seconds
    [weighted_freq_structure] = calculate_interaction_frequency(subsampled_composition_matrix,subsampled_ZOI_matrix,subsampled_composition_table.Subject,weight_option);
    [weighted_sims_structure] = permute_composition(subsampled_composition_table,subsampled_ZOI_matrix,num_sims,weight_option,replace_option,'family','subject');
    [fighandle,deltaIF_structure, pvals] = plot_simulation_results(weighted_freq_structure,weighted_sims_structure,num_sims,fignum);
end


%% Compare targets of antagonism
% Figure 2B,C,D
shuffle_option = 'both'; % shuffle_option: rows, cols, or both to shuffle antagonizers and/or baits

%Phylogroup
fignum = 0;
if fignum
    % [fighandle p fishertable] = plot_fishers_exact(interactions_by_lineage_sepi,'Phylogroup','Phylogroup',fignum);
    [simulation_structure] = permute_interactions_by_group(interactions_by_lineage_sepi,'Phylogroup',shuffle_option,1000,fignum);
    [fighandle] = plot_between_group_differences(simulation_structure,'Phylogroup','Phylogroup',fignum);
end

%Agr Type
fignum = 0;
if fignum
    % [fighandle p fishertable] = plot_fishers_exact(interactions_by_lineage_sepi,'Agr_Type','agr Type',fignum);
    [simulation_structure] = permute_interactions_by_group(interactions_by_lineage_sepi,'Agr_Type',shuffle_option,1000,fignum);
    [fighandle] = plot_between_group_differences(simulation_structure,'Agr_Type','agr Type',fignum);
end
%
%Species
fignum = 0;
if fignum
    % [fighandle p fishertable] = plot_fishers_exact(interactions_by_lineage,'Species','Species',fignum);
    [simulation_structure] = permute_interactions_by_group(interactions_by_lineage,'Species',shuffle_option,1000,fignum);
    [fighandle] = plot_between_group_differences(simulation_structure,'Species','Species',fignum);
end

%% Make heatmap to compare two families
% Figure 3C

fignum=0;
if fignum
    % create duplicated ZOI data for the following subjects
    samples_to_include = {'1AA1','1AA3','1AA4','1PA1','1PA2','1PA4','1PB2','1PB3','5PA1','5PB3'};
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
    interactions_with_replication = subsample_interaction_structure(interactions_by_lineage,included_idxs);
    for n=1:num
        interactions_with_replication.metadata(n).Name = included_labels{n};
        interactions_with_replication.metadata(n).Subject = replicated_subject_labels{n};
    end
    fighandle = overlay_interaction_heatmap(interactions_with_replication,'Name',fignum,'Subject');
end


%% Make heatmap to compare all families
fignum=0;
if fignum
    % create duplicated ZOI data for the following subjects
    samples_to_include = {'1AA1','1AA3','1AA4','1PA1','1PA2','1PA4','1PB2','1PB3','2AA1','2AA3','2PA1','2PA3','2PB1','4AA1','4AB1','5PA1','5PB3','7AA1','7AA4','7AB4','7PA1','8AA4','8AB4','8AC4','8PA3','8PB1','8PB4'};
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
    interactions_with_replication = subsample_interaction_structure(interactions_by_lineage,included_idxs);
    for n=1:num
        interactions_with_replication.metadata(n).Name = included_labels{n};
        interactions_with_replication.metadata(n).Subject = replicated_subject_labels{n};
    end
    fighandle = overlay_interaction_heatmap(interactions_with_replication,'Name',fignum,'Subject');
end

%% Are antagonizer lineages at higher abundance?
% Figure 3D

fignum = 0; % makes 2 figs
if fignum
    % remove composition data for unrepresented lineages
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
    [abundance_structure] = calculate_abundance_correlation(subsampled_composition_matrix,subsampled_ZOI_matrix);
    [fighandle] = plot_abundance_vs_antagonism(abundance_structure,'antagonism',fignum);
    [fighandle] = plot_abundance_vs_antagonism(abundance_structure,'sensitivity',fignum+1);
end

%% Make heatmap to show idiosyncratic isolates
% Figure 4A

fignum=0;
if fignum
    % create duplicated ZOI data for the following lineages
    lineages_to_include = [7,12,21,35,57,20,37,58];
    idiosyncratic_idxs = [];
    replicated_subject_labels = {};
    num = 0;
    for l = 1:numel(lineages_to_include)
        idxs = find([interactions_by_isolate_sepi.metadata.Lineage] == lineages_to_include(l));
        idiosyncratic_idxs = [idiosyncratic_idxs idxs];

    end
    interactions_idiosyncratic = subsample_interaction_structure(interactions_by_isolate_sepi,idiosyncratic_idxs);
    fighandle = overlay_interaction_heatmap(interactions_idiosyncratic,'Name',fignum,'Lineage');
end


%% Plot fold enrichment of different gene categories
% Supplementary Figure 1
fignum=0; %makes 2 figs
if fignum
    accessory_stats = readtable("figure_S1_accessory_stats_output.csv");
    gainloss_stats = readtable("figure_S1_gainloss_stats_output.csv");
    fighandle = plot_core_fold_enrichment(fignum,accessory_stats);
    fighandle = plot_core_fold_enrichment(fignum+1,gainloss_stats);
end

%% Cluster interactions by similarity and plot dendrograms
% Supplementary Figure 6
fignum=0; % clustergram is buggy, opens multiple figs
if fignum
    fighandle = plot_interaction_clustergram(interactions_by_lineage,fignum);
end

%% Mechanism upset plot
% Supplementary Figure 7
fignum=0;
if fignum
    fighandle = plot_mechanism_upset('mechanism_upset_data.csv',fignum);
end

%% Plot composition by agr type and phylogroup
% Supplementary Figure 9
fignum=0;
if fignum
    fighandle = plot_agr_and_phylo_composition(interactions_by_lineage_sepi,rep_lineages,composition_table,composition_matrix,fignum);
end

%% Plot MIC for relavant isolates
% Figure 4B
fignum = 0;
if fignum
    fighandle = plot_MIC_change('figure_4_MICs_with_revertant.xlsx',fignum);
end

% Supplementary Figure 11
fignum=0;
if fignum
    fighandle = summarize_MIC('figure_S11_MICs_summary.xlsx',fignum);
end

%% Plot growth rate per lineage
% Supplementary Figure 12B
% Do antagonizers have different growth rates?
fignum = 0; %creates 3 figs
if fignum
    figure(fignum)
    clf(fignum)
    growth_filename = "figure_S12_growth_rates_2308_2312";
    [fighandle] = plot_growthrate_vs_antagonism(interactions_by_lineage_sepi,interactions_by_isolate,composition_table,growth_filename,fignum);
end

%% Plot growth rate per lineage
% Supplementary Figure 12B
% Do antagonizers have different growth rates?
fignum = 0; %creates 3 figs
if fignum
    figure(fignum)
    clf(fignum)
    growth_filename = "figure_S12_growth_rates_2308_2312";
    [fighandle] = plot_growthrate_vs_antagonism(interactions_by_lineage_sepi,interactions_by_isolate,composition_table,growth_filename,fignum);
end

%% Plot growth rate vs abundance
fignum = 0; %creates 3 figs
if fignum
    figure(fignum)
    clf(fignum)
    growth_filename = "figure_S12_growth_rates_2308_2312";
    [fighandle] = plot_growthrate_vs_abundance(interactions_by_lineage_sepi,interactions_by_isolate,composition_table,growth_filename,fignum);
end

%% Plot growth rate of idiosyncratic isolates
% Supplementary Figure 12D
fignum = 0;
if fignum
    growth_filename = "figure_S12_idiosyncratic_growth";
    fighandle = plot_idiosyncratic_growth(growth_filename,fignum);
end

%% Are antagonizer lineages more likely to be shared?
% Supplementary Figure 12C
fignum = 0; 
if fignum
    figure(fignum)
    clf(fignum)
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
    fighandle = plot_sharing_vs_antagonism(subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix,fignum);
end




%% Other helpful analyses
list_option = false;

% Permute "same timepoint" label for each subject
fignum = 0;
if fignum
    plot_per_sample_interaction_frequency_difference(interactions_by_lineage_sepi,composition_table,groupby_option,shuffle_option,num_sims,fignum)
end

% Look for cases of potential resistance
fignum=0;
if fignum
    fighandle = plot_intralineage_variation(interactions_by_isolate,500,list_option,fignum);
end

% Do antagonizers affect shannon diversity?
fignum = 0;
if fignum
    comp_idxs = 1:size(composition_table,1);
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_by_lineage_sepi, comp_idxs);
    fighandle = plot_shannon_vs_antagonism(subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix,fignum);
end

% Do children have higer AF?
fignum=0;
if fignum
    freq = weighted_freq_structure.per_subject_interaction_freq;
    AF_child = freq(~contains(weighted_freq_structure.subjects,'P'));
    AF_adult = freq(contains(weighted_freq_structure.subjects,'P'));
    [H,P] = ttest2(AF_child,AF_adult)
end

% Compare dMRCA for antagonists vs non-antagonists
fignum = 0;
if fignum
    fighandle = plot_dMRCA_vs_antagonism(interactions_by_lineage_sepi,'dMRCA_table.csv');
end

do_stats=1;
list_option=1;
if do_stats
    % Ratio of shared / unshared lineages
    sharing_table = grpstats(composition_matrix,composition_table{:,3},"mean");
    num_shared = sum(sum(sharing_table>0)>1);
    num_unshared = sum(sum(sharing_table>0)>0)-num_shared;
    if list_option
        disp([num2str(num_shared,3) ' lineages are shared across family members and ' num2str(num_unshared,3) ' lineages are not.'])
    end
    % Superantagonizers vs nonanatagonizers
    num_ants = sum(ZOI_matrix);
    is_superant = num_ants>(numel(num_ants)./2);
    superant_fraction = sum(num_ants(is_superant))./sum(num_ants);
    nonant_fraction = sum(num_ants==0)./numel(num_ants);
    if list_option
        disp([num2str(sum(is_superant),3) ' superantagonist lineages are responsible for ' num2str(superant_fraction,3) ' of antagonisms.'])
        disp([num2str(nonant_fraction,3) ' of lineages do not antagonize any others.']) 
        disp([num2str(sum(num_ants)./(numel(ZOI_matrix)-numel(num_ants)),3) ' of interactions are antagonistic.'])
        disp(['Weighted Antagonism Frequency across the cohort is ' num2str(weighted_freq_structure.expected_freq,3)])
    end    
end


