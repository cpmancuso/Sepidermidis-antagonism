function plot_background_intensity(images,params,default_spot_locations,manually_corrected_spot_results,automated_background_results)

% Unpack variables

all_intense_meds = automated_background_results.all_intense_meds;
all_intense_means = automated_background_results.all_intense_means;
all_intense_stds = automated_background_results.all_intense_stds;
all_intense_serrs = automated_background_results.all_intense_serrs;
all_background_circle = automated_background_results.all_background_circle;
all_background_square = automated_background_results.all_background_square;

num_images = numel(images);



% Step A: produce histogram to compare all background measurements in one
% image. e.g. spot-to-spot variation per image
per_image_figure = figure(1);
clf(1)
for i=1:num_images
    subplot(8,12,i)
    histogram(all_background_circle(i,:),BinEdges=0:15:255)
    hold on
    histogram(all_background_square(i,:),BinEdges=0:15:255)
    xlim([0 255])
end

per_image_heatmap = figure(2);
clf(2)
mean_background = nanmean(all_background_circle');
heatmap(mean_background);

% Step B: produce histogram to compare all background measurements for one
% spot. e.g. image-to-image variation per spot
per_spot_figure = figure(3);
clf(3)
for n=1:96 % Assumes 96-well plate
    subplot(8,12,n)
    histogram(all_background_circle(1:num_images,n),BinEdges=0:15:255)
    hold on
    histogram(all_background_square(1:num_images,n),BinEdges=0:15:255)
    xlim([0 255])
end

per_spot_heatmap = figure(4);
clf(4)
mean_background = reshape(nanmean(all_background_circle),12,8)';
heatmap(mean_background);