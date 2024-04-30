function [fighandle, deltaIF_structure, pvals] = plot_simulation_results(freq_structure,simulation_structure,num_sims,fignum)

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [fighandle.Position(1) fighandle.Position(2) 300 300];

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



% calculate deltaInteractionFrequency
% This value measures whether antagonistic interactions are enriched or
% depleted relative to the expected interaction frequency.

per_subject_deltaIF = per_subject_interaction_freq - expected_freq;
per_subject_deltaIF_sims = per_subject_interaction_freq_sims - expected_freq_sims;

num_subjects = size(per_subject_interaction_freq_sims,2);
for s=1:num_subjects
    subplot(2,num_subjects,s)

    % histogram(per_subject_interaction_freq_sims(:,s),0:0.01:1,'Orientation','Horizontal')
    cdfplot(per_subject_deltaIF_sims(:,s))
    hold on

    permutation_val(s,1) = sum(per_subject_deltaIF(s)>=per_subject_deltaIF_sims(:,s))./num_sims;
    plot(per_subject_deltaIF(s),permutation_val(s,1),'.','MarkerSize',mksz)
    title(subjects{s})
end


families = {'1','2','4','5','7','8'};
for f=1:numel(families)
    per_family_deltaIF(f) = mean(per_subject_deltaIF(contains(subjects,families{f})));
    per_family_deltaIF_sims(:,f) = mean(per_subject_deltaIF_sims(:,contains(subjects,families{f})),2);
end

mean_per_subject_deltaIF_sims = mean(per_subject_deltaIF_sims,2);
mean_per_subject_deltaIF = mean(per_subject_deltaIF);
subject_perm_val = (1+sum(mean_per_subject_deltaIF_sims<=mean_per_subject_deltaIF))./num_sims;


mean_per_family_deltaIF_sims = mean(per_family_deltaIF_sims,2);
mean_per_family_deltaIF = mean(per_family_deltaIF);
family_perm_val = (1+sum(mean_per_family_deltaIF_sims<=mean_per_family_deltaIF))./num_sims;

subplot(2,num_subjects,[s+1:s+s])
histogram(mean_per_subject_deltaIF_sims,'FaceColor','b')
hold on
plot(mean_per_subject_deltaIF,0,'o','MarkerFaceColor','k','MarkerEdgeColor','k')
% histogram(mean_per_family_deltaIF_sims,'FaceColor','r')
% plot(mean_per_family_deltaIF,0,'o','MarkerFaceColor','b','MarkerEdgeColor','k')
title(['Mean by subject p=' num2str(subject_perm_val)]) % '; Mean by family p=' num2str(family_perm_val)])
xlabel('Mean Delta Interaction Freq')

deltaIF_structure.mean_per_subject_deltaIF_sims = mean_per_subject_deltaIF_sims;
deltaIF_structure.mean_per_subject_deltaIF = mean_per_subject_deltaIF;
deltaIF_structure.mean_per_family_deltaIF_sims = mean_per_family_deltaIF_sims;
deltaIF_structure.mean_per_family_deltaIF = mean_per_family_deltaIF;
pvals = [subject_perm_val,family_perm_val];