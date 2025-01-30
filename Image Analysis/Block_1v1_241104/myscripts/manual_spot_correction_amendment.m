function output_filename = manual_spot_correction_amendment(output_filename,images,params,default_spot_locations,manually_corrected_spot_results,list_of_images_to_correct)

% Unpack structures
all_x_coords = default_spot_locations.all_x_coords;
all_y_coords = default_spot_locations.all_y_coords;

% Correct spot locations if necessary for these images
all_x_coords = all_x_coords + 0;
all_y_coords = all_y_coords - 150;

crop_spot_dectect = params.crop_spot_dectect;
object_range = params.object_range;

all_centers = manually_corrected_spot_results.all_centers;
all_radii = manually_corrected_spot_results.all_radii;
spot_metrics = manually_corrected_spot_results.spot_metrics;
all_adj_centers = manually_corrected_spot_results.all_adj_centers;

num_images = numel(list_of_images_to_correct);

%% Step A: NA

%% Step B: Manually correct spots
for j=1:num_images
    i=list_of_images_to_correct(j); %
    figure(1)
    clf(1)
    imshow(images{i})
    title(['Manually correcting spots that could not be automatically detected in image ' num2str(i) '...'])
    drawnow

    for n=1:96 % 96-well plate
        plot_centers(n,[1,2]) = all_adj_centers(i,n,[1 2]);
        plot_radii(n)=max([all_radii(i,n),3]); % minimum radius to plot is 3, for visualization
    end
    % plot circles around each spot so far
    figure(1)
    viscircles(gca,plot_centers,plot_radii,'Color','r','LineWidth',0.5)
    hold on
    % Collect points that need to be corrected
    button=1;
    bad_spots=false(1,96);    
    title(['Image ' num2str(i) ': Click any spots to correct manually then press enter to continue...'])
    while button<=3
        [x,y,button] = ginput(1);
        if x
            [~,idx] = min(abs(x-all_x_coords)+abs(y-all_y_coords)) ;
            plot(all_x_coords(idx),all_y_coords(idx),'o','MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
            bad_spots(idx)=true;
        end
    end

    % Correct spots manually
    for n=1:96
        if bad_spots(n)
            crop_img = imcrop(images{i},[all_x_coords(n)-crop_spot_dectect/2,all_y_coords(n)-crop_spot_dectect/2,crop_spot_dectect,crop_spot_dectect]);
            [center,radius,metric,] = userfind_spot(crop_img,crop_spot_dectect,object_range);
            all_centers{i,n} = center;
            all_radii(i,n) = radius;
            spot_metrics(i,n)= metric;
        end
    end

    % Save adjusted centers individual spots
    for n=1:96
        %save to 3d matrix
        all_adj_centers(i,n,[1 2]) = all_centers{i,n}+[all_x_coords(n) all_y_coords(n)]-[crop_spot_dectect crop_spot_dectect]./2;
    end
end

title('Saving manually corrected spot results...')
disp('Saving manually corrected spot results...')
pause(1)

% Convert to structure
manually_corrected_spot_results.all_centers = all_centers;
manually_corrected_spot_results.all_radii = all_radii;
manually_corrected_spot_results.spot_metrics = spot_metrics;
manually_corrected_spot_results.all_adj_centers = all_adj_centers;

save(output_filename,"manually_corrected_spot_results")
close all