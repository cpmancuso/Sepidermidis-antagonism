function output_filename = automated_background_detection(output_filename,images,image_names,params,default_spot_locations,manually_corrected_spot_results)

% Unpack structures
object_range = params.object_range;
spot_mask_radius = params.spot_mask_radius;
rect_mask_radius = params.rect_mask_radius;
crop_radius = params.crop_radius;
int_radius = params.int_radius;

all_radii = manually_corrected_spot_results.all_radii;
all_adj_centers = manually_corrected_spot_results.all_adj_centers;

num_images = numel(images);

%% Step A: Mask spots and plate edge
% This step masks spots and plate edges and saves to file.

% Initialize
all_intense_meds = cell(num_images,96);
all_intense_means = cell(num_images,96);
all_intense_stds = cell(num_images,96);
all_intense_serrs = cell(num_images,96);
all_background_circle = nan(num_images,96);
all_background_square = nan(num_images,96);

expected_time = int8(num_images/10); %rough estimate of 6 seconds per image
background_detection_figure = figure(1);


for i=1:num_images
    tic;

    imsize = size(images{i});
    disp(['Masking Image ' num2str(i) ' of ' num2str(num_images) '. ' num2str(expected_time) ' minutes remaining...'])

    rect_masked_img = images{i};

    % Mask plate edges
    max_col=ceil(max(all_adj_centers(i,:,1))+rect_mask_radius);
    min_col=ceil(min(all_adj_centers(i,:,1))-rect_mask_radius);
    max_row=ceil(max(all_adj_centers(i,:,2))+rect_mask_radius);
    min_row=ceil(min(all_adj_centers(i,:,2))-rect_mask_radius);
    rect_mask = true(size(rect_masked_img));
    rect_mask([min_row:max_row],[min_col:max_col]) = false;
    rect_masked_img(rect_mask)=NaN; %sets to 0 since int8 doesn't take NaN    

    imwrite(rect_masked_img,['Masked Images\' image_names(i).name(1:end-4) '_rect_masked.jpg'],'jpg')

    % Mask all spots, save mask locations to add back in
    spot_masked_img = rect_masked_img;
    for n=1:96        
        [xx,yy] = ndgrid((1:imsize(1))-all_adj_centers(i,n,2),(1:imsize(2))-all_adj_centers(i,n,1));
        mask{n} = (xx.^2 + yy.^2)<spot_mask_radius^2; %logical        
        spot_masked_img(mask{n}) = NaN;
    end
    imwrite(spot_masked_img,['Masked Images\' image_names(i).name(1:end-4) '_spot_masked.jpg'],'jpg')
    imshow(spot_masked_img)

    disp(['Measuring Image ' num2str(i) ' of ' num2str(num_images) '. ' num2str(expected_time) ' minutes remaining...'])
    drawnow
    title(['Measuring Image ' num2str(i) ' of ' num2str(num_images) '. ' num2str(expected_time) ' minutes remaining...'])

%% Step B: Measure intensity
% This step measures the intensity along the radial distance away from the 
% spot edge. 

    parfor n=1:96
        %Crop images to measure intensity and background

        % Image containing spots and ZOI region, but plate edges masked
        cropped_spot_img = imcrop(rect_masked_img,[all_adj_centers(i,n,1)-crop_radius,all_adj_centers(i,n,2)-crop_radius,2*crop_radius,2*crop_radius]);
        % Image with plate edges, spots, and ZOI regions masked
        cropped_back_img = imcrop(spot_masked_img,[all_adj_centers(i,n,1)-crop_radius,all_adj_centers(i,n,2)-crop_radius,2*crop_radius,2*crop_radius]);

        % Measure intensity
        [intense_med,intense_mean,intense_std,intense_serr] = calc_radialintensity(cropped_spot_img,int_radius);
        [background_circle,background_square] = calc_background(cropped_back_img,crop_radius);
        all_intense_meds{i,n}=intense_med;
        all_intense_means{i,n}=intense_mean;
        all_intense_stds{i,n}=intense_std;
        all_intense_serrs{i,n}=intense_serr;
        all_background_circle(i,n)=background_circle;
        all_background_square(i,n)=background_square;
    end
    elapsed_time = toc;
    expected_time = round(mean([expected_time,((num_images-i)*elapsed_time./60)]));
end

%% Save outputs
title('Saving automated background results...')
disp('Saving automated background results...')
pause(1)

% Convert to structure
automated_background_results.all_radii = all_radii;
automated_background_results.all_adj_centers = all_adj_centers;
automated_background_results.all_intense_meds = all_intense_meds;
automated_background_results.all_intense_means = all_intense_means;
automated_background_results.all_intense_stds = all_intense_stds;
automated_background_results.all_intense_serrs = all_intense_serrs;
automated_background_results.all_background_circle = all_background_circle;
automated_background_results.all_background_square = all_background_square;

save(output_filename,"automated_background_results")
close(background_detection_figure)
