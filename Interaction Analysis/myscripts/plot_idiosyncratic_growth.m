function fighandle = plot_idiosyncratic_growth(growth_filename,fignum);

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [100 100 350 350];

% Options
plot_option = 'median'; % 'median' or 'allpoints'

growth_data = readtable(growth_filename);

isolates = [215,230,211, 101, 105, 313, 3131, 3133, 227, 311, 353]';

for i=1:numel(isolates)
    idxs = find(growth_data.Colony==isolates(i));
    % normalize
    gr{i} = mean([growth_data.GrowthRate(idxs(1:4)),growth_data.GrowthRate(idxs(5:8))]');
    cap{i} = mean([growth_data.SaturationOD(idxs(1:4)),growth_data.SaturationOD(idxs(5:8))]');
    yyaxis left
    set(gca,'YColor','k');
    scatter(i.*[1 1 1 1]-0.2,gr{i},'k')
    hold on
    yyaxis right
    set(gca,'YColor','k');
    scatter(i.*[1 1 1 1]+0.2,cap{i},'r')
end
xticks(1:numel(isolates))
xticklabels({'20.1','20.2','20.3','37.1','37.2','37.3','37.3-r1','37.3-r2','58.1','58.2','58.3'})
xlim([0 numel(isolates)+1])

xlabel('Isolate')
yyaxis left
ylim([0 1.2])
ylabel('Growth Rate')
yyaxis right
ylim([0 1.2])
ylabel('Saturation OD')
