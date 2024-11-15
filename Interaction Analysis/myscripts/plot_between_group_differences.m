function [fighandle] = plot_between_group_differences(simulation_structure,groupby_option,label,fignum)
fighandle = figure(fignum);
clf(fignum)
fighandle.Position=[fighandle.Position(1) fighandle.Position(2) 1200 400]; %approx 2in tall
formatSpec = '%.3f';

%% Plot Fisher's Exact
subplot(1,3,1)

same_int = nnz(simulation_structure.same_group_interaction_vector);
diff_int = nnz(simulation_structure.diff_group_interaction_vector);
same_nonint = numel(simulation_structure.same_group_interaction_vector)-nnz(simulation_structure.same_group_interaction_vector);
diff_nonint = numel(simulation_structure.diff_group_interaction_vector)-nnz(simulation_structure.diff_group_interaction_vector);

fishertable = table([same_int;diff_int],[same_nonint;diff_nonint],'RowNames',{['Same ' label],['Diff ' label]},'VariableNames',{'Antagonism','Neutral'});
[h,p,stats] = fishertest(fishertable);

% Convert to relative abundance
rel_fish_data = fishertable{:,:}'./sum(fishertable{:,:}');
rel_fish_n = sum(fishertable{:,:}');

% Plot Bar
% b=bar(rel_fish_data','Stacked');


% Plot Line
l=plot(rel_fish_data(1,:),'-o','MarkerFaceColor','w','MarkerEdgeColor','k','Color','k','LineWidth',2);
xlim([0.5 2.5])
ylim([0 0.25])

ylabel('Frequency of Antagonism')
xticks([1 2])
xticklabels({['Same (' num2str(rel_fish_n(1)) ')'],['Different (' num2str(rel_fish_n(2)) ')']})
xtickangle(0)
xlabel(label)
title(['Fisher''s Exact p = ' num2str(p,formatSpec)])
pbaspect([1 1 1])


%% Plot permutation test results
subplot(1,3,2)

observed_deltaIF = rel_fish_data(1,2) - rel_fish_data(1,1);
simulated_deltaIF = simulation_structure.diff_group_interaction_freq_sims - simulation_structure.same_group_interaction_freq_sims;
histogram(simulated_deltaIF)
hold on
scatter(observed_deltaIF,0,'r','filled')
xlim([-0.1 0.1])
pbaspect([1 1 1])
xlabel('Different Group - Same Group')

pval_high = sum(observed_deltaIF>simulated_deltaIF)./numel(simulated_deltaIF);
pval_low = sum(observed_deltaIF<simulated_deltaIF)./numel(simulated_deltaIF);
pval = min(pval_high,pval_low)*2;
title(['Permutation 2-sided p-value: ' num2str(pval)])

%% Plot AUC comparison for nonzero data
subplot(1,3,3)

same_group_interaction_size_vector = nonzeros(simulation_structure.same_group_interaction_vector);
diff_group_interaction_size_vector = nonzeros(simulation_structure.diff_group_interaction_vector);

swarmchart(ones(size(same_group_interaction_size_vector)),same_group_interaction_size_vector)
hold on
swarmchart(2*ones(size(diff_group_interaction_size_vector)),diff_group_interaction_size_vector)
plot([1 2],[median(same_group_interaction_size_vector) median(diff_group_interaction_size_vector)])

legend({['Same ' label],['Diff ' label]})
[hks,pks,statsks] = kstest2(same_group_interaction_size_vector,diff_group_interaction_size_vector);
title(['KS-test p = ' num2str(pks,formatSpec)])
ylabel('AUC of ZOI')
xticks([1 2])
xticklabels({['Same (' num2str(rel_fish_n(1)) ')'],['Different (' num2str(rel_fish_n(2)) ')']})
xlabel(label)
pbaspect([1 1 1])

