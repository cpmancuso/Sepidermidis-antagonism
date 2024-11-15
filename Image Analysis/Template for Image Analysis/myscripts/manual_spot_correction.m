function output_filename = manual_spot_correction(output_filename,images,params,default_spot_locations,automated_spot_results)


% Unpack structures
all_x_coords = default_spot_locations.all_x_coords;
all_y_coords = default_spot_locations.all_y_coords;

crop_spot_dectect = params.crop_spot_dectect;
object_range = params.object_range;

all_centers = automated_spot_results.all_centers;
all_radii = automated_spot_results.all_radii;
spot_metrics = automated_spot_results.spot_metrics;

num_images = numel(images);

%% Step A: Identify spots to be systematically excluded
% This step can be used to reduce time spent on blanks, do not use to
% indicate failure to grow on a single plate.
figure(1)
clf(1)
imshow(images{1})
% Collect blank locations to be excluded from analysis
button=1;
spots_to_exclude=false(1,96);

title('First, click any locations to systematically exclude blanks from all plates, then press enter to continue...')
drawnow
hold on

while button<=3
    [x,y,button] = ginput(1);
    if x
        [~,idx] = min(abs(x-all_x_coords)+abs(y-all_y_coords)) ;
        plot(all_x_coords(idx),all_y_coords(idx),'o','MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
        spots_to_exclude(idx)=true;
    end
end

% Populate 

%% Step B: Manually correct spots
for i=1:num_images
    figure(1)
    clf(1)
    imshow(images{i})
    title(['Manually correcting spots that could not be automatically detected in image ' num2str(i) '...'])
    drawnow

    for n=1:96 % 96-well plate
        if spots_to_exclude(n)
            all_centers{i,n} = [round(crop_spot_dectect/2) round(crop_spot_dectect/2)];
            all_radii(i,n) = 0;
            spot_metrics(i,n)= 0;
            continue
        elseif spot_metrics(i,n)<0.05 %if no center or obviously bad center was found
            crop_img = imcrop(images{i},[all_x_coords(n)-crop_spot_dectect/2,all_y_coords(n)-crop_spot_dectect/2,crop_spot_dectect,crop_spot_dectect]);
            [center,radius,metric] = userfind_spot(crop_img,crop_spot_dectect,object_range);
            all_centers{i,n} = center;
            all_radii(i,n) = radius;
            spot_metrics(i,n)= metric;
        end
        % save adjusted center relative to image for plotting
        plot_centers(n,[1,2]) = all_centers{i,n}+[all_x_coords(n) all_y_coords(n)]-[crop_spot_dectect crop_spot_dectect]./2;
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