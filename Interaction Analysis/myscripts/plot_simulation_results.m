function [fighandle, deltaAF_structure, pvals] = plot_simulation_results(freq_structure,simulation_structure,num_sims,fignum)

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [fighandle.Position(1) fighandle.Position(2) 200 200];

plot_CDFs = false;

% unpack structures
expected_freq = freq_structure.expected_freq;
per_sample_interaction_freq = freq_structure.per_sample_interaction_freq;
per_subject_interaction_freq = freq_structure.per_subject_interaction_freq;
subjects = freq_structure.subjects;
expected_freq_sims = simulation_structure.expected_freq_sims;
per_sample_interaction_freq_sims = simulation_structure.per_sample_interaction_freq_sims;
per_subject_interaction_freq_sims = simulation_structure.per_subject_interaction_freq_sims;

% check sizes
if numel(subjects)~=size(per_subject_interaction_freq_sims,2) | size(per_subject_interaction_freq,2)~=size(per_subject_interaction_freq_sims,2)
    error('Number of subjects in data and simulations must match.')
end

mksz = 20;

% calculate deltaAF
% This value measures whether antagonistic interactions are enriched or
% depleted relative to the expected interaction frequency.

per_subject_deltaAF = per_subject_interaction_freq - expected_freq;
per_subject_deltaAF_sims = per_subject_interaction_freq_sims - expected_freq_sims;

num_subjects = size(per_subject_interaction_freq_sims,2);

%% Optionally plot CDFs for each subject
if plot_CDFs
    for s=1:num_subjects
        subplot(2,num_subjects,s)
        cdfplot(per_subject_deltaAF_sims(:,s))
        hold on
    
        permutation_val(s,1) = sum(per_subject_deltaAF(s)>=per_subject_deltaAF_sims(:,s))./num_sims;
        plot(per_subject_deltaAF(s),permutation_val(s,1),'.','MarkerSize',mksz)
        title(subjects{s})
    end
    % prep for permutation plot results
    subplot(2,num_subjects,[s+1:s+s])
end

% calculate per_family deltaAF
families = {'1','2','4','5','7','8'};
for f=1:numel(families)
    per_family_deltaAF(f) = mean(per_subject_deltaAF(contains(subjects,families{f})));
    per_family_deltaAF_sims(:,f) = mean(per_subject_deltaAF_sims(:,contains(subjects,families{f})),2);
end

mean_per_subject_deltaAF_sims = mean(per_subject_deltaAF_sims,2);
mean_per_subject_deltaAF = mean(per_subject_deltaAF);
subject_perm_val = (1+sum(mean_per_subject_deltaAF_sims<=mean_per_subject_deltaAF))./num_sims;

mean_per_family_deltaAF_sims = mean(per_family_deltaAF_sims,2);
mean_per_family_deltaAF = mean(per_family_deltaAF);
family_perm_val = (1+sum(mean_per_family_deltaAF_sims<=mean_per_family_deltaAF))./num_sims;

% Plot deltaAF simulations and observed
histogram(mean_per_subject_deltaAF_sims,'FaceColor',[0.7 0.7 0.7])
hold on
plot(mean_per_subject_deltaAF,0,'o','MarkerFaceColor','k','MarkerEdgeColor','k')
title(['Mean by subject p=' num2str(subject_perm_val)])
xlabel('Mean Delta Antagonism Freq')
ylabel('Simulations')
pbaspect([2 1 1])

% Save resulting data structure
deltaAF_structure.mean_per_subject_deltaAF_sims = mean_per_subject_deltaAF_sims;
deltaAF_structure.mean_per_subject_deltaAF = mean_per_subject_deltaAF;
deltaAF_structure.mean_per_family_deltaAF_sims = mean_per_family_deltaAF_sims;
deltaAF_structure.mean_per_family_deltaAF = mean_per_family_deltaAF;
pvals = [subject_perm_val,family_perm_val];