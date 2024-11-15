function output_filename = automated_spot_detection(output_filename,images,params,default_spot_locations)

% Unpack variables
all_x_coords = default_spot_locations.all_x_coords;
all_y_coords = default_spot_locations.all_y_coords;

crop_spot_dectect = params.crop_spot_dectect;
object_range = params.object_range;

num_images = numel(images);
expected_time = int8(num_images/3); %rough estimate of 20 seconds per image

spot_detection_figure = figure(1);
for i=1:numel(images)
    tic
    % Display image
    imshow(images{i})
    title(['Analysing Image ' num2str(i) ' of ' num2str(num_images) '. ' num2str(expected_time) ' minutes remaining...'])
    drawnow
    img = images{i};
    % Find center and radius for each spot
    parfor n=1:96 % assumes 96-well plate, needs to match 
        disp(['Analysing Image ' num2str(i) ' spot ' num2str(n)])

        crop_img = imcrop(img,[all_x_coords(n)-crop_spot_dectect/2,all_y_coords(n)-crop_spot_dectect/2,crop_spot_dectect,crop_spot_dectect]);

        [center,radius,metric] = imfind_spot(crop_img,crop_spot_dectect,object_range);
        all_centers{i,n} = center;
        all_radii(i,n) = radius;
        spot_metrics(i,n)= metric;
    end
    elapsed_time = toc;
    expected_time = round(mean([expected_time,((num_images-i)*elapsed_time./60)]));
end

title('Saving automated spot results...')
disp('Saving automated spot results...')
pause(1)

% Convert to structure
automated_spot_results.all_centers = all_centers;
automated_spot_results.all_radii = all_radii;
automated_spot_results.spot_metrics = spot_metrics;

save(output_filename,"automated_spot_results")
close(spot_detection_figure)
