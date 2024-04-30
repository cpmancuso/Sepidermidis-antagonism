function [simulation_structure] = permute_composition(composition_table,ZOI_matrix,num_sims,weight_option,replace_option,pool_option,shuffle_option)

% weight_option: weighted or nonweighted, for calculating interaction frequency
% replace_option: replace or noreplace
% pool_option: family or population, for determining range of permutation
% shuffle_option: column or subject, for determining method of permutation


composition_matrix = composition_table{:,4:end};
num_samples = size(composition_table,1);
unique_subjects = unique(composition_table.Subject);
num_subjects = numel(unique_subjects);
num_lineages = size(composition_matrix,2);


% initialize
expected_freq_sims = zeros(num_sims,1);
per_sample_interaction_freq_sims = zeros(num_sims,num_samples);
per_subject_interaction_freq_sims = zeros(num_sims,num_subjects);
if ~islogical(replace_option)
    replace_bool = strcmp(replace_option,'replace');
end

% permute based on permutation option, first be pool option, then by
% shuffle_option

switch pool_option
    case 'population'
        switch shuffle_option
            case 'column'
                for n=1:num_sims
                    sim_composition_matrix = composition_matrix;
                    perm_idxs = datasample(1:num_lineages,num_lineages,'Replace',replace_bool);
                    sim_composition_matrix = composition_matrix(:,perm_idxs);
                    sim_freq_structure = calculate_interaction_frequency(sim_composition_matrix,ZOI_matrix,composition_table.Subject,weight_option);
                    expected_freq_sims(n) = sim_freq_structure.expected_freq;
                    per_sample_interaction_freq_sims(n,:) = sim_freq_structure.per_sample_interaction_freq;
                    per_subject_interaction_freq_sims(n,:) = sim_freq_structure.per_subject_interaction_freq;
                end
            case 'subject'
                for n=1:num_sims
                    sim_composition_matrix = composition_matrix;
                    for s=1:num_subjects
                        perm_idxs = datasample(1:num_lineages,num_lineages,'Replace',replace_bool);
                        idxs = strcmp({composition_table.Subject{:}},unique_subjects{s});
                        sim_composition_matrix(idxs,:) = composition_matrix(idxs,perm_idxs);
                    end
                    sim_freq_structure = calculate_interaction_frequency(sim_composition_matrix,ZOI_matrix,composition_table.Subject,weight_option);
                    expected_freq_sims(n) = sim_freq_structure.expected_freq;
                    per_sample_interaction_freq_sims(n,:) = sim_freq_structure.per_sample_interaction_freq;
                    per_subject_interaction_freq_sims(n,:) = sim_freq_structure.per_subject_interaction_freq;
                end
        end
    case 'family'
        families = {'1','2','4','5','7','8'};
        switch shuffle_option
            case 'column'
                for n=1:num_sims
                    sim_composition_matrix = composition_matrix;
                    for f=1:numel(families)
                        % find indexes of subjects and lineages that belong to this family
                        fam_idxs = [composition_table.Family{:}]==str2num(families{f});
                        fam_lineages = find(sum(composition_matrix(fam_idxs,:))>0); % find lineages present in composition
                        % permute over only subjects and lineages of this family  
                        perm_idxs = datasample(1:numel(fam_lineages),numel(fam_lineages),'Replace',replace_bool);
                        sim_composition_matrix(fam_idxs,fam_lineages) = composition_matrix(fam_idxs,fam_lineages(perm_idxs));
                    end
                    sim_freq_structure = calculate_interaction_frequency(sim_composition_matrix,ZOI_matrix,composition_table.Subject,weight_option);
                    expected_freq_sims(n) = sim_freq_structure.expected_freq;
                    per_sample_interaction_freq_sims(n,:) = sim_freq_structure.per_sample_interaction_freq;
                    per_subject_interaction_freq_sims(n,:) = sim_freq_structure.per_subject_interaction_freq;
                end
            case 'subject'
                for n=1:num_sims
                    sim_composition_matrix = composition_matrix;
                    for f=1:numel(families)
                        % find indexes of subjects and lineages that belong to this family
                        fam_idxs = [composition_table.Family{:}]==str2num(families{f});
                        fam_lineages = find(sum(composition_matrix(fam_idxs,:))>0); % find lineages present in composition
                        fam_subjects = unique({composition_table.Subject{fam_idxs}});
                        % permute over only subjects and lineages of this family  
                        for s=1:numel(fam_subjects)
                            perm_idxs = datasample(1:numel(fam_lineages),numel(fam_lineages),'Replace',replace_bool);
                            idxs = strcmp({composition_table.Subject{:}},fam_subjects{s});
                            sim_composition_matrix(idxs,fam_lineages) = composition_matrix(idxs,fam_lineages(perm_idxs));
                        end
                    end
                    sim_freq_structure = calculate_interaction_frequency(sim_composition_matrix,ZOI_matrix,composition_table.Subject,weight_option);
                    expected_freq_sims(n) = sim_freq_structure.expected_freq;
                    per_sample_interaction_freq_sims(n,:) = sim_freq_structure.per_sample_interaction_freq;
                    per_subject_interaction_freq_sims(n,:) = sim_freq_structure.per_subject_interaction_freq;
                end
        end
end

simulation_structure.expected_freq_sims = expected_freq_sims;
simulation_structure.per_sample_interaction_freq_sims = per_sample_interaction_freq_sims;
simulation_structure.per_subject_interaction_freq_sims = per_subject_interaction_freq_sims;