clc
clear all
close all
cwd = pwd();
addpath([pwd '/myscripts'])
load('params.mat')
imagedir = 'C:\Users\cmanc\Dropbox\Lieberman Lab\Personal lab notebooks\Chris Mancuso\Data\220712_Block2_ZOI_assay\2022_07_15';
spots_name = ['spots_' imagedir(end-9:end) '_measured_10-Aug-2022'];
data_name = ['data_' imagedir(end-9:end) '_measured_10-Aug-2022'];
processed_name = ['processed_data_' imagedir(end-9:end) '_measured_16-Nov-2022'];

old_spots = 'oldspots_2022_07_15.mat';
old_data = '';

images = dir(imagedir);
images=images(contains({images.name},'JPG'));
imgprefix = images(1).name(1:end-6);
mkdir('Masked Images')

% %% Load images and calibrate to image of choice
% if isfile('calibration.mat')
%     load('calibration.mat')
% else
%     i = 1; %image to calibrate with
%     img = imread_orient([images(i).folder '\' images(i).name]);
%     imshow(img)
%     set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
%     corners = calibrate_platemap(gcf);
%     calx = sort(corners(1,:)); %sort to correct for clicking order
%     caly = sort(corners(2,:)); %sort to correct for clicking order
%     xvals = linspace(mean(calx([1,2])),mean(calx([3,4])),12);
%     yvals = linspace(mean(caly([1,2])),mean(caly([3,4])),8);
%     xall = repmat(xvals,1,8);
%     yall = repelem(yvals,12);
%     clear calx caly corners img
%     save('calibration.mat','xall','yall','xvals','yvals')
% end
% 
% 
%load images
disp('Loading Images...')
parfor i=1:numel(images)
    img{i} = imread_orient([images(i).folder '\' images(i).name]);
end
% 
% %parameters
% crop_dim = 350; %diameter of image to look for spot in
% object_range = [60 120]; %size range of spots
% expected_time = 30;
% spot_mask_radius = 180; %size of spots to mask 
% rect_mask_radius = 150; %edge around outer rows of spots to keep
% crop_radius = 300; %radius of image to measure intensity on
% int_radius = 200; %radius to measure intensity on
% 
% %% Loop through each image to find centers and edges
% %process images, comment block if using pre-identified spots
% for i=1:numel(images)
%     tic
%     %show image
%     figure(1)
%     imshow(img{i})
%     title(['Analysing Image ' num2str(i) ' of 96. ' num2str(expected_time) ' minutes remaining...'])
%     drawnow
%     parfor n=1:96 %find center and radius from cropped image
%         disp(['Analysing Image ' num2str(i) ' spot ' num2str(n)])
%         
%         crop_img{n} = imcrop(img{i},[xall(n)-crop_dim/2,yall(n)-crop_dim/2,crop_dim,crop_dim]);
% 
%         [center,radius,metric] = imfind_spot(crop_img{n},object_range);
%         centers{i,n} = center;
%         all_radii(i,n) = radius;
%         metrics(i,n)= metric;
%     end
%     elapsed_time = toc;
%     expected_time = round(mean([expected_time,((96-i)*elapsed_time./60)]));
% end
% 
% %% User input section, initial image processing
% for i=1:numel(images)
% 
%     figure(1)
%     imshow(img{i})
%     title(['Analysing Image ' num2str(i) '...'])
%     drawnow
%     
%     for n=1:96
%         if metrics(i,n)<0.05 %if no center or obviously bad center was found
%             crop_img{n} = imcrop(img{i},[xall(n)-crop_dim/2,yall(n)-crop_dim/2,crop_dim,crop_dim]);
%             [center,radius,metric] = userfind_spot(crop_img{n},object_range);
%             centers{i,n} = center;
%             all_radii(i,n) = radius;
%             metrics(i,n)= metric;
%         end
%         % save adjusted center relative to image for plotting
%         plot_centers(n,[1,2]) = centers{i,n}+[xall(n) yall(n)]-[crop_dim crop_dim]./2;
%         plot_radii(n)=max([all_radii(i,n),3]);
%     end
%     % plot circles for each ZOI so far
%     figure(1)
%     viscircles(gca,plot_centers,plot_radii,'Color','r','LineWidth',0.5)
%     hold on
%     % Collect points that need to be corrected
%     button=1;
%     bad_spots=false(1,96);    
%     title(['Image ' num2str(i) ': Click any spots to correct manually then press enter to continue...'])
%     while button<=3
%         [x,y,button] = ginput(1);
%         if x
%             [~,idx] = min(abs(x-xall)+abs(y-yall)) ;
%             plot(xall(idx),yall(idx),'o','MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
%             bad_spots(idx)=true;
%         end
%     end
%     
%     % Correct spots manually
%     for n=1:96
%         if bad_spots(n)
%             crop_img{n} = imcrop(img{i},[xall(n)-crop_dim/2,yall(n)-crop_dim/2,crop_dim,crop_dim]);
%             [center,radius,metric,] = userfind_spot(crop_img{n},object_range);
%             centers{i,n} = center;
%             all_radii(i,n) = radius;
%             metrics(i,n)= metric;
%         end
%     end
%     
%     % Save adjusted centers individual spots
%     for n=1:96
%         %save to 3d matrix
%         all_adj_centers(i,n,[1 2]) = centers{i,n}+[xall(n) yall(n)]-[crop_dim crop_dim]./2;
%     end
% end
% 
% % Save data
% disp('Saving data...')
% save([spots_name '.mat'],'all_radii','all_adj_centers');
% disp(['Saved to ' spots_name '.mat'])
% 
% 
% %% Alternatively, load old data structure
% % process images, comment block if identifying de-novo spots
% load(old_spots)
% if exist('all_centers','var')==1 %old data saved adjusted center in cells
%     all_adj_centers = zeros(96,96,2);
%     old_radii = all_radii;
%     all_radii = zeros(96,96);
%     for i=1:96
%         for n=1:96
%             all_adj_centers(i,n,[1 2]) = all_centers{i}(n,:);
%             all_radii(i,n) = old_radii{i}(n);
%         end
%     end
% end
% save([spots_name '.mat'],'all_radii','all_adj_centers');    

% %% Fix centers if necessary
% load([spots_name '.mat'])
% load('calibration.mat')
% for i=89 %enter image to correct directly
%     close all
%     figure(1)
%     imshow(img{i})
%     hold on
%     title(['Analysing Image ' num2str(i) '...'])
%     plot_centers = reshape(all_adj_centers(i,:,:),96,2);
%     plot_radii = all_radii(i,:);
%     viscircles(gca,plot_centers,plot_radii,'Color','r','LineWidth',0.5)
%     plot_centers = reshape(all_adj_centers(i,:,:),96,2);
%     plot_radii = spot_mask_radius.*ones(96,1);
%     viscircles(gca,plot_centers,plot_radii,'Color','b','LineWidth',0.5)
%     drawnow
%     % Collect points that need to be corrected
%     button=1;
%     bad_spots=false(1,96);    
%     title(['Image ' num2str(i) ': Click any spots to correct manually then press enter to continue...'])
%     while button<=3
%         [x,y,button] = ginput(1);
%         if x
%             [~,idx] = min(abs(x-xall)+abs(y-yall)) ;
%             plot(xall(idx),yall(idx),'o','MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
%             bad_spots(idx)=true;
%         end
%     end
%     % Correct spots manually
%     for n=1:96
%         if bad_spots(n)
%             crop_img{n} = imcrop(img{i},[xall(n)-crop_dim/2,yall(n)-crop_dim/2,crop_dim,crop_dim]);
%             [center,radius,metric,] = userfind_spot(crop_img{n},object_range);
%             all_radii(i,n) = radius;
%             all_adj_centers(i,n,[1 2]) = center+[xall(n) yall(n)]-[crop_dim crop_dim]./2;
%         end
%     end
%     close all
% end
% save([spots_name '.mat'],'all_radii','all_adj_centers');    
% 
% %% Using finalized centers, mask spots to measure background
% correcting = 0;
% if correcting
%     load([data_name '.mat'])
% else
%     cropped_spot_img = cell(96,96);
%     cropped_back_img = cell(96,96);
%     all_intense_meds = cell(96,96);
%     all_intense_means = cell(96,96);
%     all_intense_stds = cell(96,96);
%     all_intense_serrs = cell(96,96);
%     all_background_circle = zeros(96,96);
%     all_background_circle = zeros(96,96);
% end
% load([spots_name '.mat'])
% close all
% cropped_spot_img = cell(96,96);
% cropped_back_img = cell(96,96);
% all_intense_meds = cell(96,96);
% all_intense_means = cell(96,96);
% all_intense_stds = cell(96,96);
% all_intense_serrs = cell(96,96);
% all_background_circle = zeros(96,96);
% all_background_circle = zeros(96,96);
% expected_time = 10;
% for i=1:96
%     tic;
% 
%     img{i} = imread_orient([images(i).folder '\' images(i).name]);
%     imsize = size(img{i});
%     disp(['Masking Image ' num2str(i) ' of 96. ' num2str(expected_time) ' minutes remaining...'])
%     
%     rect_masked_img = img{i};
%     
%     % Mask plate edges
%     max_col=ceil(max(all_adj_centers(i,:,1))+rect_mask_radius);
%     min_col=ceil(min(all_adj_centers(i,:,1))-rect_mask_radius);
%     max_row=ceil(max(all_adj_centers(i,:,2))+rect_mask_radius);
%     min_row=ceil(min(all_adj_centers(i,:,2))-rect_mask_radius);
%     rect_mask = true(size(rect_masked_img));
%     rect_mask([min_row:max_row],[min_col:max_col]) = false;
%     rect_masked_img(rect_mask)=NaN; %sets to 0 since int8 doesn't take NaN    
% 
%     imwrite(rect_masked_img,['Masked Images\' images(i).name(1:end-4) '_rect_masked.jpg'],'jpg')
%     
%     % Mask all spots, save mask locations to add back in
%     spot_masked_img = rect_masked_img;
%     for n=1:96        
%         [xx,yy] = ndgrid((1:imsize(1))-all_adj_centers(i,n,2),(1:imsize(2))-all_adj_centers(i,n,1));
%         mask{n} = (xx.^2 + yy.^2)<spot_mask_radius^2; %logical        
%         spot_masked_img(mask{n}) = NaN;
%     end
%     imwrite(spot_masked_img,['Masked Images\' images(i).name(1:end-4) '_spot_masked.jpg'],'jpg')
%     imshow(spot_masked_img)
%     title(['Measuring Image ' num2str(i)])
%     drawnow
%     disp(['Measuring Image ' num2str(i)])
%     parfor n=1:96
%         %Crop images to measure intensity and background
%         cropped_spot_img = imcrop(rect_masked_img,[all_adj_centers(i,n,1)-crop_radius,all_adj_centers(i,n,2)-crop_radius,2*crop_radius,2*crop_radius]);
%         cropped_back_img = imcrop(spot_masked_img,[all_adj_centers(i,n,1)-crop_radius,all_adj_centers(i,n,2)-crop_radius,2*crop_radius,2*crop_radius]);
% 
%         % Measure intensity
%         [intense_med,intense_mean,intense_std,intense_serr] = calc_radialintensity(cropped_spot_img,int_radius);
%         [background_circle,background_square] = calc_background(cropped_back_img,crop_radius);
%         all_intense_meds{i,n}=intense_med;
%         all_intense_means{i,n}=intense_mean;
%         all_intense_stds{i,n}=intense_std;
%         all_intense_serrs{i,n}=intense_serr;
%         all_background_circle(i,n)=background_circle;
%         all_background_square(i,n)=background_square;
%     end
%     elapsed_time = toc;
%     expected_time = round(mean([expected_time,((96-i)*elapsed_time./60)]));
% end
% 
% % Save data
% disp('Saving data...')
% save([data_name '.mat'],'all_radii','all_adj_centers','all_intense_meds','all_intense_means','all_intense_stds','all_intense_serrs','all_background_circle','all_background_square');
% disp(['Saved to ' data_name '.mat'])
% save('params.mat','imgprefix','crop_dim', 'object_range', 'spot_mask_radius', 'rect_mask_radius','crop_radius','int_radius')

%% Process data for plotting
clearvars -except spots_name data_name processed_name imgprefix
load([data_name '.mat'])
load('params.mat')
thresh = 8;

for i=1:96
    for n=1:96
        % calculate normalized intensities
        norm_x(i,n,1:int_radius)=(1:int_radius)-all_radii(i,n)+1; %shift relative to radius
        norm_int1(i,n,1:int_radius)=all_intense_meds{i,n}-all_background_circle(i,n); %smaller background region seems more accurate, generally lower
        norm_int2(i,n,1:int_radius)=all_intense_meds{i,n}-all_background_square(i,n);
        if norm_int1(i,n,int_radius)<-thresh
            norm_int3(i,n,1:int_radius)=norm_int1(i,n,1:int_radius); %large ZOI, can't local correct
        else
            norm_int3(i,n,1:int_radius)=norm_int1(i,n,1:int_radius)-mean(norm_int1(i,n,int_radius-4:int_radius)); %take edge value as background
        end        
        
        roi = norm_x(i,n,:) > 0; %logical
        %find longest zoi roi using runlength
        zoi_roi = ZOI_runlength(norm_int3(i,n,:),roi);
        
        ZOI_width(i,n) = sum(zoi_roi,3); %pixels
        ZOI_depth(i,n) = -min(norm_int3(i,n,roi)); %intensity
        
        ZOI_area(i,n) = -trapz(norm_int3(i,n,roi)); %all periphery, invert
        ZOI_subarea(i,n) = -trapz(norm_int3(i,n,zoi_roi)); %just ZOI, invert
        if (ZOI_depth(i,n)>thresh || ZOI_subarea(i,n)>1000) && (all_radii(i,n)>10);
            ZOI_call(i,n)=ZOI_subarea(i,n);
        else
            ZOI_call(i,n)=0;
        end
    end
end

%save data
disp('Saving data...')
save([processed_name '.mat'],'all_radii','all_adj_centers','all_intense_meds','all_intense_means','all_intense_stds','all_intense_serrs','all_background_circle','all_background_square','norm_x','norm_int1','norm_int2','norm_int3','ZOI_width','ZOI_depth','ZOI_call','ZOI_area','ZOI_subarea');
disp(['Saved to ' processed_name '.mat'])
