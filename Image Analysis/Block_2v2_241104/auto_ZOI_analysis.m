clc
clear all
close all
cwd = pwd();
addpath([pwd '/myscripts'])

%% Choose which steps to complete

% Choose image directory
imagedir = 'C:\Users\cmanc\Dropbox\Lieberman Lab\Personal lab notebooks\Chris Mancuso\Data\2024_10_18_block2v2';
if isdir(imagedir)
    image_names = dir(imagedir);
    image_names=image_names(contains({image_names.name},'JPG')); % to find images later
    num_images = numel(image_names);
    imgprefix = image_names(1).name(1:end-6);
    analysis_suffix = strcat('analyzed_',date); % suffix to append to intermediate analysis files
else
    disp([imagedir ' not found...'])
end

% User can provide a filename here to bypass this step, or leave it empty
% to rerun analysis for all subsequent steps.
% example_step1 = 'example.mat'
% example_step2 = ''

% define default spot locations (same for all plates
default_spot_locations_file = 'default_spot_locations.mat';

% automate detection of spot edges 
automated_spot_results_file = 'automated_spot_results_2024_10_18_block2v2_analyzed_05-Nov-2024.mat';

% manually correct automated spot detection
manually_corrected_spot_results_file = 'manually_corrected_spot_results_2024_10_18_block2v2_analyzed_05-Nov-2024.mat';


% automate measurement of background intensity
automated_background_results_file = 'automated_background_results_2024_10_18_block2v2_analyzed_05-Nov-2024.mat';

% calculate ZOIs using background intensities
automated_ZOI_results_file = 'automated_ZOI_results_2024_10_18_block2v2_analyzed_08-Nov-2024.mat';

% Image analysis parameters, do not change between experiments, only between imaging set up
params_filename = 'params.mat';
if ~isfile(params_filename)
    params.crop_spot_dectect = 350; % diameter of image to look for spot in
    params.object_range = [60 120]; % size range of spots
    params.spot_mask_radius = 180; % size of spots to mask 
    params.rect_mask_radius = 150; % edge around outer rows of spots to keep
    params.crop_radius = 300; % radius to crop image to measure intensity on
    params.int_radius = 200; % radius to measure intensity on
    params.plate_location = [700 350 4800 3200]; % relevant area of plate, used for cropping
    params.local_background_correction_threshold = 8; % intensity difference used to find ZOIs too large to local correct  
    save(params_filename,'params')
else
    disp(['Loaded image parameters from file...'])
    load(params_filename)
end

%% Generate run report

rerun_all_subsequent = false;

run_define_default_spot_locations = ~isfile(default_spot_locations_file)||rerun_all_subsequent;
if run_define_default_spot_locations
    rerun_all_subsequent = true;
    disp('No default spot location file exists, will define default spot locations...')
else
    disp(['Loading default spot locations from: ' default_spot_locations_file])
end

run_automated_spot_detection = ~isfile(automated_spot_results_file)||rerun_all_subsequent;
if run_automated_spot_detection
    rerun_all_subsequent = true;
    disp('No automated spot results file exists, will run automated spot detection...')
else
    disp(['Loading automated spot locations from: ' automated_spot_results_file])
end

run_manual_spot_correction = ~isfile(manually_corrected_spot_results_file)||rerun_all_subsequent;
if run_manual_spot_correction
    rerun_all_subsequent = true;
    disp('No manually corrected spot results file exists, will run manual spot correction...')
else
    disp(['Loading manually corrected spot locations from: ' manually_corrected_spot_results_file])
end

run_automated_background_detection = ~isfile(automated_background_results_file)||rerun_all_subsequent;
if run_automated_background_detection
    rerun_all_subsequent = true;
    disp('No automated background results file exists, will run automated background detection...')
else
    disp(['Loading automated background results from: ' automated_background_results_file])
end

run_automated_ZOI_calculation = ~isfile(automated_ZOI_results_file)||rerun_all_subsequent;
if run_automated_ZOI_calculation
    rerun_all_subsequent = true;
    disp('No automated ZOI results file exists, will run automated ZOI detection...')
else
    disp(['Loading automated ZOI results from: ' automated_ZOI_results_file])
end

tempfolders = split(imagedir,'\');
expt_name = tempfolders{end}; %used for filenaming

%% STEP 1: Define Default Spot Locations
% This step asks the user to define the location of the four corner wells
% in order to define default spot locations for all images to be analyzed

if run_define_default_spot_locations
    output_filename = 'default_spot_locations.mat'; %to save results
    cal_img = 1; %image number to calibrate with
    default_spot_locations_file = define_default_spot_locations(image_names,cal_img,output_filename,params.plate_location);
    load(default_spot_locations_file)
else
    load(default_spot_locations_file)
end

%% Load All Images
disp(['Loading ' num2str(num_images) ' images...'])
parfor i=1:numel(image_names)
    images{i} = imread_orient([image_names(i).folder '\' image_names(i).name],params.plate_location);
end

%% STEP 2: Automated Spot Detection
% This step loads each image and automatically detects the edge of all 96
% spots on the plate. Centers, radii, and quality metrics are saved for
% each spot on each image.

if run_automated_spot_detection
    output_filename = ['automated_spot_results_' expt_name '_' analysis_suffix]; %to save results
    automated_spot_results_file = automated_spot_detection(output_filename,images,params,default_spot_locations);
    load(automated_spot_results_file)
else
    load(automated_spot_results_file)
end

%% STEP 3: Manual Spot Correction
% This step shows the user each image that has been analyzed automatically.
% Any automatically flagged spots will require user input, but additional
% spots can also be added for review. 

if run_manual_spot_correction
    output_filename = ['manually_corrected_spot_results_' expt_name '_' analysis_suffix]; %to save results
    manually_corrected_spot_results_file = manual_spot_correction(output_filename,images,params,default_spot_locations,automated_spot_results);
    load(manually_corrected_spot_results_file)
else
    load(manually_corrected_spot_results_file)
end


%% STEP 3.5: Manual Spot Correction Amendment
% This step permits the user to correct specific plate images for any spots
% which were miscalled.

list_of_images_to_correct = [];

if ~isempty(list_of_images_to_correct)
    output_filename = ['manually_corrected_spot_results_' expt_name '_' analysis_suffix]; %to save results
    manually_corrected_spot_results_file = manual_spot_correction_amendment(output_filename,images,params,default_spot_locations,manually_corrected_spot_results,list_of_images_to_correct);
    load(manually_corrected_spot_results_file)
else
    load(manually_corrected_spot_results_file)
end

%% STEP 4: Automated background detection
% This step masks spots and plate edges in order to measure the background
% around each spot. It then measures the intensity at each radial distance
% from each spot.


if run_automated_background_detection
    mkdir('Masked Images')
    output_filename = ['automated_background_results_' expt_name '_' analysis_suffix]; %to save results
    automated_background_results_file = automated_background_detection(output_filename,images,image_names,params,default_spot_locations,manually_corrected_spot_results);
    load(automated_background_results_file)
else
    load(automated_background_results_file)
end

%% STEP 5: Review background intensity variation 
% This step plots local background intensity for each spot. 
plot_background_intensity(images,params,default_spot_locations,manually_corrected_spot_results,automated_background_results)

%% STEP 6: Automated ZOI calculation
% This step calculates the width, depth, area, and subarea of the region
% near the spot where ZOIs would appear. An initial call is made for the
% ZOI, which is later filtered in the Interaction Analysis script.

if run_automated_ZOI_calculation
    output_filename = ['automated_ZOI_results_' expt_name '_' analysis_suffix]; %to save results
    automated_ZOI_results_file = automated_ZOI_calculation(output_filename,images,params,default_spot_locations,manually_corrected_spot_results,automated_background_results);
    load(automated_ZOI_results_file)
else
    load(automated_ZOI_results_file)
end