function fighandle = plot_interaction_frequency_for_each_family(composition_table,interactions_by_lineage,fignum)

fighandle = figure(fignum);
clf(fignum)


families = {'1','2','4','5','7','8'};
all_expect = {};
all_per_sample = {};
all_per_subject = {};
p = 0;
for f=1:numel(families)
    %filter composition and find lineages to track for ZOI
    fam_idxs = [composition_table.Family{:}]==str2num(families{f});
    [subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix] = subsample_composition(composition_table,interactions_by_lineage, fam_idxs);
    
    [family_freq_structure] = calculate_interaction_frequency(subsampled_composition_matrix,subsampled_ZOI_matrix,subsampled_composition_table.Subject,'weighted');
    all_expect{f} = family_freq_structure.expected_freq;
    all_per_sample{f} = family_freq_structure.per_sample_interaction_freq;
    all_per_subject{f} = family_freq_structure.per_subject_interaction_freq;
    
    subplot(2,6,f)
    num = numel(all_per_sample{f});
    xvals = (1:num);
    sh = stem(xvals,all_per_sample{f},'BaseValue',all_expect{f},'MarkerFaceColor','b','Color','k');
    xticks(xvals)
    xticklabels(subsampled_composition_table.SID(:))
    xtickangle(90)
    xlim([0 num+1])
    ylim([0 0.5])
    sh.BaseLine.Color = 'r';

    
    subplot(2,6,f+6)
    num = numel(all_per_subject{f});
    xvals = (1:num);
    sh = stem(xvals,all_per_subject{f},'BaseValue',all_expect{f},'MarkerFaceColor','b','Color','k');
    xticks(xvals)
    xticklabels(family_freq_structure.subjects)
    xtickangle(90)
    xlim([0 num+1])
    ylim([0 0.5])
    sh.BaseLine.Color = 'r';

end

subplot(2,6,1)
ylabel('Within Sample Interaction Frequency')

subplot(2,6,7)
ylabel('Within Subject Interaction Frequency')