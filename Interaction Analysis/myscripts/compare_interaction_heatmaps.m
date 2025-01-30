function fighandle = compare_interaction_heatmaps(interaction_structure,label_option,fignum,overlay_option)

% params for plotting
num = numel(interaction_structure.metadata);
labels = {interaction_structure.metadata.(label_option)};
fig_position = [10.6000 1009 2560 1.3168e+03]; %specific to monitor setup

% create figures
fighandle = figure(fignum);
fighandle.Position = fig_position;
clf(fighandle);

% calculate correlations
ZOI_jaccard = pdist([interaction_structure.ZOI_call1(:),interaction_structure.ZOI_call2(:)]','jaccard');
disp(['Jaccard similarity is: ' num2str(1-ZOI_jaccard)])

ZOI_SMC = sum(interaction_structure.ZOI_call1(:)==interaction_structure.ZOI_call2(:))./numel(interaction_structure.ZOI_call(:));
disp(['Simple matching coefficient is: ' num2str(ZOI_SMC)])

ZOI_1v2 = sum((interaction_structure.ZOI_call1(:)>0)&(interaction_structure.ZOI_call2(:)>0))./sum((interaction_structure.ZOI_call1(:)>0));
disp(['Percent of condition 1 antagonisms also observed in condition 2:' num2str(ZOI_1v2)])

ZOI_2v1 = sum((interaction_structure.ZOI_call1(:)>0)&(interaction_structure.ZOI_call2(:)>0))./sum((interaction_structure.ZOI_call2(:)>0));
disp(['Percent of condition 2 antagonisms also observed in condition 1:' num2str(ZOI_2v1)])

% Create plot data
plot_data = zeros(size(interaction_structure.ZOI_call));
plot_data(interaction_structure.ZOI_call1&~interaction_structure.ZOI_call2) = 1;
plot_data(interaction_structure.ZOI_call2&~interaction_structure.ZOI_call1) = 2;
plot_data(interaction_structure.ZOI_call1&interaction_structure.ZOI_call2) = 3;


% Define categorical colormap
hexColors = {'#FAD149', '#3870B8','#58A051'};
categoryRGB = zeros(length(hexColors), 3);
for i = 1:length(hexColors)
    hex = hexColors{i};
    categoryRGB(i, :) = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
end
categoryRGB = [0.1 0.1 0.1; categoryRGB];
colormap(categoryRGB);

% Make heatmap 
plot_data = plot_data'; %transpose for imagesc
image_handle = imagesc(plot_data,[0 3]);
colorbar
hold on
xlabel('Reciever Lawn')
ylabel('Producer Spot')
xticks(1:num)
xticklabels(labels)
h=gca; h.XAxis.TickLength = [0 0];
xtickangle(90)
yticks(1:num)
yticklabels(labels)
h=gca; h.YAxis.TickLength = [0 0];

% Make grid
xrange = [1 num];
dx = diff(xrange)/(num-1);
xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,num+1);
gridhandle = mesh(xg,xg,zeros(num+1));
gridhandle.FaceColor = 'none';
gridhandle.EdgeColor = 'k';

%% Draw overlays
if ~strcmp(overlay_option,'none')
    
    % Find group boundaries
    if isnumeric([interaction_structure.metadata.(overlay_option)])
        group_labels = unique([interaction_structure.metadata.(overlay_option)],'stable');
        boundaries = zeros(size(group_labels));
        for n=1:numel(group_labels)
            boundaries(n) = find([interaction_structure.metadata.(overlay_option)]==group_labels(n), 1);
        end
    else
        group_labels = unique(cellstr({interaction_structure.metadata.(overlay_option)}),'stable');
        boundaries = zeros(size(group_labels));
        for n=1:numel(group_labels)
            boundaries(n) = find(strcmp(cellstr({interaction_structure.metadata.(overlay_option)}),group_labels(n)), 1);
        end
    end
    boundaries = [boundaries num+1];

    % White dashes
    [horz,vert] = meshgrid(boundaries);
    grouphandle = mesh(horz-0.5,vert-0.5,zeros(size(horz)));
    grouphandle.FaceColor = 'none';
    grouphandle.EdgeColor = 'w';
    grouphandle.LineStyle = '--';

    % Red boxes
    for g=2:numel(boundaries)
        s = boundaries(g-1);
        e = boundaries(g);
        patch('Faces',[1:4],'Vertices',[s e; e e; e s; s s]-0.5,'EdgeColor','r','FaceColor','none','LineWidth',2)
    end
end
pbaspect([1 1 1])

