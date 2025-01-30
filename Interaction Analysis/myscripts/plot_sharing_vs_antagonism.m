function [fighandle,max_ants_per_sample] = plot_sharing_vs_antagonism(subsampled_composition_table,subsampled_composition_matrix,subsampled_ZOI_matrix,fignum)

axis_option = 'fraction';
% Note, this function assumes that composition matrix and ZOI matrix have
% been subsampled such that ZOI matrix only contains lineages that are
% present in composition matrix. 

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [fighandle.Position(1) fighandle.Position(2) 900 300];

num_samples = size(subsampled_composition_matrix,1);
num_lineages = size(subsampled_composition_matrix,2);
num_subjects = numel(unique({subsampled_composition_table.Subject{:}}));
mean_antagonism_frequency = mean(subsampled_ZOI_matrix);
num_ants = nansum(subsampled_ZOI_matrix);

npairs = 0;
labels = {};
for s=1:num_samples
    sample_lineages = subsampled_composition_matrix(s,:)>0;
    max_ants_per_sample(s)=max(num_ants(sample_lineages));
    % Consider each pair to see if same family, different subject
    for f=s:num_samples
        is_family_member = subsampled_composition_table.Family{s}==subsampled_composition_table.Family{f} & ~strcmp(subsampled_composition_table.Subject{s},subsampled_composition_table.Subject{f});
        if is_family_member
            family_member_lineages = subsampled_composition_matrix(f,:)>0;
            npairs=npairs+1;
            max_antagonism(npairs) = max(mean_antagonism_frequency(sample_lineages|family_member_lineages)); % this sample or partner
            num_shared_lineages(npairs) = sum(sample_lineages&family_member_lineages);
            frac_shared_lineages(npairs) = num_shared_lineages(npairs)./sum(sample_lineages|family_member_lineages,2);
            labels{npairs} = [subsampled_composition_table.SID{s} ':' subsampled_composition_table.SID{f}];
        end
    end
end

%reshape
max_antagonism = max_antagonism';
num_shared_lineages = num_shared_lineages';
frac_shared_lineages = frac_shared_lineages';


lineage_x = max_antagonism;
lineage_y = num_shared_lineages;
sample_x = max_antagonism;
sample_y = num_shared_lineages;
color_x = '#C03026'; %red
color_y = '#3870B8'; %blue

% color_x = '#C03026'; %red
% color_y = '#808285'; %grey
% color_x = '#58A051'; %green
% color_y = '#3870B8'; %blue
% color_y = '#8C67AC'; %purple


[linear_rho,linear_p] = corr(sample_x,sample_y,'Type','Pearson');
[spearman_rho,spearman_p] = corr(lineage_x,lineage_y,'Type','Spearman');

t = tiledlayout(1, 3);

%% Linear correlation, samples
nexttile;
mdl = fitlm(sample_x,sample_y);
scatter(sample_x,sample_y,20,'k','filled','MarkerFaceAlpha',0.5)
hold on
h=plot(mdl);
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
delete(dataHandle);
fitHandle = findobj(h,'DisplayName','Fit');
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);

xlim([0 1])
% ylim([0 1])
pbaspect([1 1 1])

xlabel('Max lineage antagonism proportion') 
ylabel('# of shared lineages')
title(['Correlation of lineage sharing and antagonism' newline 'Linear rho = ' num2str(linear_rho,'%.2f') ' p = ' num2str(linear_p,'%.3f')])
legend('hide')

%% Pearson correlation, bars
nexttile([1 2]);
[sorted_x,idxs] = sort(lineage_x);
sorted_y = lineage_y(idxs);
bar(-sorted_x,'FaceColor',color_x);
hold on
bar(sorted_y./max(sorted_y),'FaceColor',color_y);

xticks(1:numel(sorted_x))
xticklabels(labels(idxs))
xlim([0 numel(sorted_x)+1])
set(gca,'box','off')

ylim([-1 1])
yticks([-1 -0.5 0 0.5 1])
yticklabels([-1 -0.5 0 max(sorted_y)/2 max(sorted_y)])
pbaspect([2 1 1])

xlabel('Lineage')
title(['Correlation of lineage sharing and antagonism' newline 'Spearman rho = ' num2str(spearman_rho,'%.2f') ' p = ' num2str(spearman_p,'%.3f')])
