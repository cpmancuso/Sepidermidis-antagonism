function interaction_structure = build_interaction_structure(datafiles)

idx_range={
[1:96],[1:96];
[97:192],[1:96];
[1:96],[97:192];
[97:192],[97:192]};

%load data into big matrix
for b=1:4
    load(datafiles{b})

    full_adj_centers(idx_range{b,1},idx_range{b,2},:) = all_adj_centers;
%     full_background_circle(idx_range{b,1},idx_range{b,2})  = all_background_circle;
%     full_background_square(idx_range{b,1},idx_range{b,2}) = all_background_square;
%     full_int_means(idx_range{b,1},idx_range{b,2}) = all_intense_means;
%     full_int_meds(idx_range{b,1},idx_range{b,2}) = all_intense_meds;
%     full_int_serrs(idx_range{b,1},idx_range{b,2}) = all_intense_serrs;
%     full_int_stds(idx_range{b,1},idx_range{b,2}) = all_intense_stds;
    full_radii(idx_range{b,1},idx_range{b,2}) = all_radii;
%     full_norm_intensities_circle(idx_range{b,1},idx_range{b,2},:) = norm_int1;
%     full_norm_intensities_square(idx_range{b,1},idx_range{b,2},:) = norm_int2;
    full_norm_int(idx_range{b,1},idx_range{b,2},:) = norm_int3; % locally-corrected for large ZOIs
    full_norm_x(idx_range{b,1},idx_range{b,2},:) = norm_x;
%     full_ZOI_area(idx_range{b,1},idx_range{b,2}) = ZOI_area;
%     full_ZOI_call(idx_range{b,1},idx_range{b,2}) = ZOI_call;
    full_ZOI_depth(idx_range{b,1},idx_range{b,2}) = ZOI_depth;
    full_ZOI_subarea(idx_range{b,1},idx_range{b,2}) = ZOI_subarea;
    full_ZOI_width(idx_range{b,1},idx_range{b,2}) = ZOI_width;
end
clear('all_adj_centers','all_background_circle','all_background_square','all_intense_means','all_intense_meds','all_intense_serrs','all_radii')
clear('norm_int1','norm_int2','norm_int3','norm_x','ZOI_area','ZOI_call','ZOI_depth','ZOI_subarea','ZOI_width')

interaction_structure = struct();
interaction_structure.centers = full_adj_centers;
interaction_structure.radii = full_radii;
interaction_structure.norm_int = full_norm_int;
interaction_structure.norm_x = full_norm_x;
interaction_structure.ZOI_depth = full_ZOI_depth;
interaction_structure.ZOI_AUC = full_ZOI_subarea;
interaction_structure.ZOI_width = full_ZOI_width;

