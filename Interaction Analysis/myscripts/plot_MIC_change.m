function fighandle = plot_MIC_change(data_filename,fignum)

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [fighandle.Position(1) fighandle.Position(2) 350 250];

range = 'A1:E7'; 
[data, headers] = xlsread(data_filename, range);
isolates = headers(2:end,1);
conditions = headers(1,2:end);

% Calculate fold change
fc_data = data;
for c =1:numel(conditions)
    fc_data(:,c) = data(:,c)./data(1,c);
end
%%
colororder({'#0E536D','#1595BA','#84B6D0','#6D6E71'})
bar(log2(fc_data),'edgecolor','none')
ylabel('log_2(Normalized MIC)')
xlabel('Isolate')
xticklabels(isolates)
legend(conditions,'location','southwest')
xlim([0 numel(isolates)+1])
ylim([-5.5 1])
xtickangle(0)

% for c =1:numel(conditions)
%     stem((1:numel(isolates))+(c-1)/(2*numel(conditions)),log2(fc_data(:,c)),'MarkerSize',3,'MarkerEdgeColor','k')
%     hold on
% end 
% ylabel('log_2(Normalized MIC)')
% xlabel('Isolate')
% xticks(1:numel(isolates))
% xticklabels(isolates)
% legend(conditions)
% xlim([0 numel(isolates)+1])