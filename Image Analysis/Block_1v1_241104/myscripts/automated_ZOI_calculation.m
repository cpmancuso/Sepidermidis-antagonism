function output_filename = automated_ZOI_calculation(output_filename,images,params,default_spot_locations,manually_corrected_spot_results,automated_background_results)

% unpack structure
object_range = params.object_range;
spot_mask_radius = params.spot_mask_radius;
rect_mask_radius = params.rect_mask_radius;
crop_radius = params.crop_radius;
int_radius = params.int_radius;
local_background_correction_threshold = params.local_background_correction_threshold;

all_radii = manually_corrected_spot_results.all_radii;
all_adj_centers = manually_corrected_spot_results.all_adj_centers;

num_images = numel(images);

all_intense_meds = automated_background_results.all_intense_meds;
all_intense_means = automated_background_results.all_intense_means;
all_intense_stds = automated_background_results.all_intense_stds;
all_intense_serrs = automated_background_results.all_intense_serrs;
all_background_circle = automated_background_results.all_background_circle;
all_background_square = automated_background_results.all_background_square;


for i=1:num_images
    for n=1:96
        % calculate normalized intensities
        norm_x(i,n,1:int_radius)=(1:int_radius)-all_radii(i,n)+1; %shift relative to radius
        norm_int_circle(i,n,1:int_radius)=all_intense_meds{i,n}-all_background_circle(i,n); %smaller background region seems more accurate, generally lower
        norm_int_square(i,n,1:int_radius)=all_intense_meds{i,n}-all_background_square(i,n);
        if norm_int_circle(i,n,int_radius)<-local_background_correction_threshold
            norm_int_local(i,n,1:int_radius)=norm_int_circle(i,n,1:int_radius); %large ZOI, can't local correct
        else
            norm_int_local(i,n,1:int_radius)=norm_int_circle(i,n,1:int_radius)-mean(norm_int_circle(i,n,int_radius-4:int_radius)); %take edge value as background
        end        

        roi = norm_x(i,n,:) > 0; %logical
        %find longest zoi roi using runlength
        zoi_roi = ZOI_runlength(norm_int_local(i,n,:),roi);

        ZOI_width(i,n) = sum(zoi_roi,3); %pixels
        ZOI_depth(i,n) = -min(norm_int_local(i,n,roi)); %intensity

        ZOI_area(i,n) = -trapz(norm_int_local(i,n,roi)); %all periphery, invert
        ZOI_subarea(i,n) = -trapz(norm_int_local(i,n,zoi_roi)); %just ZOI, invert
        if (ZOI_depth(i,n)>local_background_correction_threshold || ZOI_subarea(i,n)>1000) && (all_radii(i,n)>0);
            ZOI_call(i,n)=ZOI_subarea(i,n);
        else
            ZOI_call(i,n)=0;
        end
    end
end

%% Save Outputs
title('Saving automated ZOI results...')
disp('Saving automated ZOI results...')
pause(1)

% Convert to structure
automated_ZOI_results.all_radii = all_radii;
automated_ZOI_results.all_adj_centers = all_adj_centers;
automated_ZOI_results.all_intense_meds = all_intense_meds;
automated_ZOI_results.all_intense_means = all_intense_means;
automated_ZOI_results.all_intense_stds = all_intense_stds;
automated_ZOI_results.all_intense_serrs = all_intense_serrs;
automated_ZOI_results.all_background_circle = all_background_circle;
automated_ZOI_results.all_background_square = all_background_square;

automated_ZOI_results.norm_x = norm_x;
automated_ZOI_results.norm_int_circle = norm_int_circle;
automated_ZOI_results.norm_int_square = norm_int_square;
automated_ZOI_results.norm_int_local = norm_int_local;

automated_ZOI_results.ZOI_area = ZOI_area;
automated_ZOI_results.ZOI_depth = ZOI_depth;
automated_ZOI_results.ZOI_width = ZOI_width;
automated_ZOI_results.ZOI_subarea = ZOI_subarea;
automated_ZOI_results.ZOI_call = ZOI_call;

save(output_filename,"automated_ZOI_results")

