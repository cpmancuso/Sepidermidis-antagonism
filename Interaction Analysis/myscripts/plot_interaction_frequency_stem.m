function [fighandle p] = plot_interaction_frequency_stem(expected_freq,per_group_interaction_freq,group_labels,groupname,fignum)
fighandle = figure(fignum);
clf(fignum)
fighandle.Position=[fighandle.Position(1) fighandle.Position(2) 300 200]; %approx 2in sq


% parameters
sort_plot = 1;

if sort_plot
    [plot_data,plot_order] = sort(per_group_interaction_freq);
else
    plot_data = per_group_interaction_freq;
    plot_order = 1:numel(per_group_interaction_freq);
end

sh = stem(plot_data,'BaseValue',expected_freq,'MarkerFaceColor','w','Color','k');
ylim([0 0.5])
ylabel('Within Sample Interaction Frequency')
xlim([0 numel(per_group_interaction_freq)+1])
xticks(1:numel(per_group_interaction_freq))
xticklabels(group_labels(plot_order))
xtickangle(90)
yticks([0 0.1 0.2 0.3 0.4 0.5])
sh.BaseLine.Color = 'r';
sh.BaseLine.LineWidth = 2;
xlabel(groupname)
set(gca,'box','off')

delta_interaction_freq = per_group_interaction_freq - expected_freq;
[h,p,ci,stats] = ttest(delta_interaction_freq);
pbaspect([2 1 1])
