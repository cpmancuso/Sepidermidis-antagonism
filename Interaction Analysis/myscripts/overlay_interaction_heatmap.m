function fighandle = overlay_interaction_heatmap(plot_structure,label_option,fignum,overlay_option)
% VARIABLES:
% plot_structure : the interaction structure to be plotted, this may be sorted, subsampled
% reference_structure : the unfiltered, original interaction structure used to point to images
% labels : metadata column chosen for labelling
% image_dir : the directories where images can be found
% fignum : figure number
% ZOI_depth_threshold : optional paramter, to show ZOI depth used for call

%FUNCTIONS:
% mouse_click(src,eventData) : respond to mouse click with plot
% KeyPressCb(src,eventData) : respond to arrow keys with plot
% display_plot(x,y) : retrieve image data and ZOI data and plot

% params for plotting
num = numel(plot_structure.metadata);
labels = {plot_structure.metadata.(label_option)};
fig_position = [10.6000 1009 2560 1.3168e+03]; %specific to monitor setup

% create figures
fighandle = figure(fignum);
fighandle.Position = fig_position;
clf(fighandle);

default_color=[1 1 1];
intensity_map=[0.1 0.1 0.1; 0.5 0.5 0.5; 0.7 0.7 0.7; 0.9 0.9 0.9];
steps = 4; %1;
colormap(intensity_map.*default_color)



% Make heatmap 
plot_data = log10(plot_structure.ZOI_AUC.*plot_structure.ZOI_call);
plot_data = plot_data'; %transpose for imagesc

image_handle = imagesc(plot_data,[0 steps]);
set(image_handle, 'AlphaData', ~isnan(plot_data))
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
set(h,'Color','#3C2415')

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
    if isnumeric([plot_structure.metadata.(overlay_option)])
        group_labels = unique([plot_structure.metadata.(overlay_option)],'stable');
        boundaries = zeros(size(group_labels));
        for n=1:numel(group_labels)
            boundaries(n) = find([plot_structure.metadata.(overlay_option)]==group_labels(n), 1);
        end
    else
        group_labels = unique(cellstr({plot_structure.metadata.(overlay_option)}),'stable');
        boundaries = zeros(size(group_labels));
        for n=1:numel(group_labels)
            boundaries(n) = find(strcmp(cellstr({plot_structure.metadata.(overlay_option)}),group_labels(n)), 1);
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

