function plot_per_sample_interaction_frequency_difference(interactions_by_lineage,composition_table,groupby_option,shuffle_option,num_sims,fignum)
fighandle = figure(fignum);
clf(fignum)
fighandle.Position=[fighandle.Position(1) fighandle.Position(2) 200 200]; %approx 2in sq

% groupby_option: metadata column to group by
% shuffle_option: rows, cols, or both to shuffle antagonizers and/or baits


formatSpec = '%.3f';
num_lineages = numel(interactions_by_lineage.metadata);
metadata_table = struct2table(interactions_by_lineage.metadata);
rep_lineages = [interactions_by_lineage.metadata.Lineage];
composition_matrix = composition_table{:,rep_lineages+3};
num_reps = size(composition_matrix,2);

% 
all_same = [];
all_diff = [];
for s=1:27
    same_group_interaction_vector = [];
    diff_group_interaction_vector = [];
    in_sample = composition_matrix(s,:)>0;
    for m=1:num_reps
        for n=1:num_reps
            if m~=n
                if in_sample(m)&&in_sample(n)
                    same_group_interaction_vector = [same_group_interaction_vector; interactions_by_lineage.ZOI_call(m,n).*interactions_by_lineage.ZOI_AUC(m,n)];
                elseif in_sample(n) && ~in_sample(m)
                    diff_group_interaction_vector = [diff_group_interaction_vector; interactions_by_lineage.ZOI_call(m,n).*interactions_by_lineage.ZOI_AUC(m,n)];
                end
            end
        end
    end
    if isempty(same_group_interaction_vector)
        same_group_interaction_vector = [0];
    end


    all_same(s) = nnz(same_group_interaction_vector)./numel(same_group_interaction_vector);
    all_diff(s) = nnz(diff_group_interaction_vector)./numel(diff_group_interaction_vector);


    plot(s,all_same(s),'ro')
    hold on
    plot(s,all_diff(s),'bo')
end
xticks(1:27)
xticklabels(composition_table.SID)

ylabel('Nonweighted interaction frequency')
legend('Sample lineages vs. coresident lineages','Sample lineages vs. nonresident lineages')


stat_vector = all_same - all_diff
[H,P] = ttest(stat_vector)


%group by sample
[group, id] = findgroups(composition_table.Subject);

for g=1:max(group)
    subject_same(g) = mean(all_same(group==g));
    subject_diff(g) = mean(all_diff(group==g));
end


[subject_diff,i] = sort(subject_diff,'ascend');
subject_same = subject_same(i);
id = id(i);


figure(fignum+1)
clf(fignum+1)
plot(subject_same,'ro')
hold on
plot(subject_diff,'bo')
xticks(1:max(group))
xticklabels(id)

ylabel('Nonweighted interaction frequency')
legend('Sample lineages vs. coresident lineages','Sample lineages vs. nonresident lineages')



