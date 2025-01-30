function output_filename = define_default_spot_locations(images,cal_img,output_filename,plate_location)

calibration_figure = figure();
img = imread_orient([images(cal_img).folder '\' images(cal_img).name],plate_location);
imshow(img)
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
title('Carefully click the center of the 4 corner wells')

x_corner_coords = [];
y_corner_coords = [];

for num = 1:4
    [x,y] = ginput(1);
    x_corner_coords(num) = x;
    y_corner_coords(num) = y;
    hold on
    plot(x,y,'o','MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
end

title('Saving inputs...')
pause(1)

% Sort to correct for clicking order
x_corner_coords = sort(x_corner_coords); %sort to correct for clicking order
y_corner_coords = sort(y_corner_coords); %sort to correct for clicking order

% Generate linearly spaced x and y coordinates
x_coord_values = linspace(mean(x_corner_coords([1,2])),mean(x_corner_coords([3,4])),12);
y_coord_values = linspace(mean(y_corner_coords([1,2])),mean(y_corner_coords([3,4])),8);
all_x_coords = repmat(x_coord_values,1,8);
all_y_coords = repelem(y_coord_values,12);

% Convert to structure
default_spot_locations.all_x_coords = all_x_coords;
default_spot_locations.all_y_coords = all_y_coords;
default_spot_locations.x_coord_values = x_coord_values;
default_spot_locations.y_coord_values = y_coord_values;

save(output_filename,"default_spot_locations")
close(calibration_figure)
