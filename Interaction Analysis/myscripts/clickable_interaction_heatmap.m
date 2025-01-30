function fighandle = clickable_interaction_heatmap(plot_structure,reference_structure,label_option,image_dir,fignum,ZOI_depth_threshold)
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
subfig_position = [227.4000 364.2000 810.4000 344]; %specific to monitor setup
input = load('params.mat','params');
crop_dim=input.params.crop_spot_dectect;
crop_radius=input.params.crop_radius;
int_radius=input.params.int_radius;

% create figures
subfignum = 123;
subfighandle = figure(subfignum); 
subfighandle.Position = subfig_position;
clf(subfighandle);

fighandle = figure(fignum);
fighandle.Position = fig_position;
clf(fighandle);

% QC
if ~issorted([reference_structure.metadata.Index])
    error('Reference interaction structure must remain sorted by index...')
elseif numel(reference_structure.metadata)~=192
    error('Reference interaction structure must remain unfiltered...')
end


default_color=[1 1 1];
steps = 4;
intensity_map=linspace(0.1,1,steps);
colormap(default_color.*intensity_map')

% Make heatmap 
plot_data = log10(plot_structure.ZOI_AUC.*plot_structure.ZOI_call);
plot_data = plot_data'; %transpose for imagesc
image_handle = imagesc(plot_data,[0 steps]);
colorbar()

image_handle.ButtonDownFcn = @mouse_click; %clickable
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
gridhandle.ButtonDownFcn = @mouse_click; %clickable

% Functions
    function mouse_click(src,eventData)
        % Get coordinates of click 
        coords = eventData.IntersectionPoint;
        global x %corresponds to reciever lawn
        global y %corresponds to producer spot
        x=round(coords(1));
        y=round(coords(2));
        display_plot(x,y)
        set(gcf,'KeyPressFcn',@KeyPressCb);
    end

    function KeyPressCb(src,eventData)
%         fprintf('key pressed: %s\n',eventData.Key);
        global x %corresponds to reciever lawn
        global y %corresponds to producer spot
        if strcmp(eventData.Key,'rightarrow')
        x=min(x+1,num);
        elseif strcmp(eventData.Key, 'leftarrow')
        x=max(x-1,1);
        elseif strcmp(eventData.Key,'uparrow') %reversed due to plot
        y=max(y-1,1);
        elseif strcmp(eventData.Key,'downarrow') %reversed due to plot
        y=min(y+1,num);
        end
        display_plot(x,y)
    end

    function display_plot(x,y)
        % Convert plot index to original index for 
        i=plot_structure.metadata(x).Index;
        n=plot_structure.metadata(y).Index;
        ZOI_call = plot_structure.ZOI_call(x,y); %NOTE: call based on group
        center = [reference_structure.centers(i,n,:)];
        radius = [reference_structure.radii(i,n)];
        norm_x = reshape([reference_structure.norm_x(i,n,:)],[],int_radius);
        norm_int = reshape([reference_structure.norm_int(i,n,:)],[],int_radius);
        ZOI_AUC = reference_structure.ZOI_AUC(i,n); %NOTE: AUC based on img
        
        % Make rectangle on selected part of heatmap
        global pt
        if exist('pt','var')
            delete(pt)
        end
        pt = rectangle('Position',[x-0.5 y-0.5 1 1],'EdgeColor','r','LineWidth',1);
        
        % Identify image set to use
        blocks = [1 3; 2 4];
        blockset = blocks(reference_structure.metadata(i).Block,reference_structure.metadata(n).Block);
        imgprefix = image_dir{blockset};
        imgnum = i;
        if i>96
            imgnum=i-96;
        else
            imgnum=i;
        end
        
        % Load image from masked image saved in folder
        figure(subfignum)
        subplot(1,2,1)
        filename = [imgprefix,sprintf('%02d',imgnum),'_rect_masked.jpg'];

        subimage = imcrop(imread(filename),[center(1)-crop_radius,center(2)-crop_radius,2*crop_radius,2*crop_radius]);
        imshow(subimage)
        hold on
        % viscircles(gca,[crop_radius+1,crop_radius+1],radius,'Color','r','LineWidth',0.1);
        
        % Plot ZOI intensity data from this image
        subplot(1,2,2)
        hold off
        plot(norm_x,norm_int,'k')
        hold on
        lw = 0.5; xr = [-100,100]; yr = [-50,75];
        plot(xr,[0,0],'r-','LineWidth',lw)
        plot(xr,[-ZOI_depth_threshold,-ZOI_depth_threshold],'r--','LineWidth',lw)
        plot([0,0],yr,'r-','LineWidth',lw)
        xlim(xr); ylim(yr)
        
        % Construct labels
        labeli = [' [' char(reference_structure.metadata(i).Species) ' ' char(num2str(reference_structure.metadata(i).Lineage)) '] '];
        labeln = [' [' char(reference_structure.metadata(n).Species) ' ' char(num2str(reference_structure.metadata(n).Lineage)) '] '];
        
%         title({['Spot ' char(reference_structure.metadata(n).Name) labeln],['Lawn ' char(reference_structure.metadata(i).Name) labeli]});
        title({['Spot ' num2str(reference_structure.metadata(n).Stock) labeln],['Lawn ' num2str(reference_structure.metadata(i).Stock) labeli]});
        legend('Intensity','','Threshold')
        xlabel('Pixels from Colony Edge')
        ylabel('Normalized Intensity')
        if ZOI_call
            text(30,50,{['ZOI AUC: ' num2str(round(ZOI_AUC))],['ZOI Call: TRUE']},'FontSize',8)
        else
            text(30,50,{['ZOI AUC: ' num2str(round(ZOI_AUC))],['ZOI Call: FALSE']},'FontSize',8)
        end
        
        % Return to main fig
        figure(fighandle)
    end
end



