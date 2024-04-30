function fighandle = plot_core_fold_enrichment(fignum,datatable)

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [100 100 700 350];

% Remove rare annotations and non-annotations
total_core = datatable{1,2};
total_noncore = datatable{1,3};
datatable = datatable(2:end,:);
total = datatable{:,2} + datatable{:,3};
datatable = datatable(total>3,:);
datatable(strcmp(datatable.Group,'Function unknown'),:)=[];
datatable(strcmp(datatable.Group,'Other COG categories'),:)=[];
datatable(strcmp(datatable.Group,'Unassigned'),:)=[];

% Sort by Odds ratio, rescale
raw_datatable = datatable;
datatable.OddsRatio = log2(datatable.OddsRatio) ;
num_entries = size(datatable,1);
for n=1:num_entries
    if datatable.OddsRatio(n)==Inf
        datatable.OddsRatio(n) = 10;
    end
end
datatable = sortrows(datatable,'OddsRatio','ascend');
% Move database values to the top
dbidxs = contains(datatable.Group, 'Database');
dbrows = datatable(dbidxs, :);
datatable(dbidxs,:) = [];
datatable = [datatable; dbrows];

barh(1:num_entries,datatable.OddsRatio)
yticks(1:num_entries)
yticklabels(datatable.Group)
xlim([-4 4])
ylim([0,num_entries+1])
xlabel('log2(Odds Ratio)')
title([datatable.Properties.VariableDescriptions{3} ' n=' num2str(total_noncore) ' vs. Core n=' num2str(total_core)])

yyaxis right
total = datatable{:,2} + datatable{:,3};
isnaive = datatable.p_value < 0.05;
plot(-3.75*ones(sum(isnaive),1),find(isnaive),'r*')
hold on
issig = datatable.p_value < 0.05./num_entries;
plot(-3.6*ones(sum(issig),1),find(issig),'r*')

yticks(1:num_entries)
yticklabels(num2str(total))
ylim([0 num_entries+1])