function fighandle = plot_mechanism_upset(data_filename,fignum)

barcolor = 0.5*ones(1,3);
fighandle = figure(fignum);
clf(fignum)
fighandle.Position=[100 100 350 350];

data = readtable(data_filename);
data = fillmissing(data,'constant',0);
data = data(:,[1:11]);
numvars = size(data,2);

% combo barplot
subplot(2,3,[2 3])
num = 0;


% expt groups
[expt_groups,~,expt_idxs] = unique(data{:,2:6},'rows');

% BGC groups
[BGC_groups,~,BGC_idxs] = unique(data{:,7:11},'rows');


for m=1:size(expt_groups,1)
    for n=1:size(BGC_groups,1)
        idxs = (expt_idxs==m)&(BGC_idxs==n);
        if sum(idxs)==0
            continue
        else
            num = num+1;
        end
        bar(num,sum(idxs),'FaceColor',barcolor)
        hold on
        if m~=n
            labels{num} = [data.Properties.VariableNames{m} '+' data.Properties.VariableNames{n}];
        else
            labels{num} = [data.Properties.VariableNames{m}];
        end
    end
end
xticks(1:num)
xlim([0 num+1])
xticklabels({})
ylabel('Count')
ylim([0 5])
yticks([0:5])
% xticklabels(labels)

% label scatter plot
subplot(2,3,[5 6])
num = 0;
for m=1:size(expt_groups,1)
    for n=1:size(BGC_groups,1)
        idxs = (expt_idxs==m)&(BGC_idxs==n);
        comb_idxs = find([expt_groups(m,:) BGC_groups(n,:)]);
        if sum(idxs)==0
            continue
        else
            num = num+1;
        end
        plot(num*ones(size(comb_idxs)),comb_idxs,'LineWidth',1,'Marker','o','Color','k','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','r')
        hold on
        if m~=n
            labels{num} = [data.Properties.VariableNames{m} '+' data.Properties.VariableNames{n}];
        else
            labels{num} = [data.Properties.VariableNames{m}];
        end
    end
end
xticks(1:num)
xlim([0 num+1])
xticklabels({})
yticklabels({})
ylim([0,numvars])

% xticklabels(labels)
% yticks(1:numvars-1)
% yticklabels({data.Properties.VariableNames{2:numvars}})
grid
set(gca, 'YDir','reverse')


%% category labels
subplot(2,3,4)
for m=2:numvars
    barh(m-1,sum(data{:,m}),'FaceColor',barcolor)
    hold on
end
ylim([0,numvars])
xlim([0 15])
set(gca, 'YDir','reverse')
set(gca, 'XDir','reverse')
yticks(1:numvars)
yticklabels({data.Properties.VariableNames{2:numvars}})
set(gca,'TickLabelInterpreter','none')
xlabel('Count')





