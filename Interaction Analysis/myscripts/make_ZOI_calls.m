function ZOI_calls = make_ZOI_calls(interaction_structure,ZOI_depth_threshold,ZOI_noise_threshold,ZOI_AUC_lower_threshold,ZOI_AUC_upper_threshold)

% Make determination call whether ZOI measurements indicate a true ZOI,
% returns logical array, can be used to mask interaction_structure
% temporarily stores the criterion for each decision
% 0: uncharacterized
% 1: legitimate ZOI, ZOI_depth_threshold criterion
% 2: legitimate ZOI, ZOI_AUC_upper_threshold criterion
% 3: legitimate ZOI, both criteria
% -1: spot too small or blank
% -2: ZOI below both qualifications, ZOI_depth_threshold & ZOI_AUC_upper_threshold
% -3: ZOI below minimum size, ZOI_AUC_lower_threshold criterion
% -4: ZOI too noisy, ZOI_noise_threshold criterion

noise_range = 150:200; %range in pixels over which to determine noise, ideally outsize expected ZOI region
[sz,~] = size(interaction_structure.ZOI_AUC); %expect square
ZOI_calls = false(sz,sz);
ZOI_criteria = int8(zeros(sz,sz));

for x = 1:sz
    for y =1:sz
        if interaction_structure.radii(x,y) < 10 %if no spot, then no ZOI
            ZOI_calls(x,y) = false;
            ZOI_criteria(x,y) = -1;
        else % is not a blank spot
            if interaction_structure.ZOI_depth(x,y)<ZOI_depth_threshold && interaction_structure.ZOI_AUC(x,y)<ZOI_AUC_upper_threshold %if fail both qualifiers
                ZOI_calls(x,y) = false;
                ZOI_criteria(x,y) = -2;
            else % has at least one qualification
                if interaction_structure.ZOI_AUC(x,y)<ZOI_AUC_lower_threshold %if below minimum size
                    ZOI_calls(x,y) = false;
                    ZOI_criteria(x,y) = -3;
                else % big enough
                    % calculate noise
                    noise = range(gradient(reshape(interaction_structure.norm_int(x,y,noise_range),1,[])));
                    if interaction_structure.ZOI_depth(x,y)<noise*ZOI_noise_threshold % too noisy
                        ZOI_calls(x,y) = false;
                        ZOI_criteria(x,y) = -4;
                    else % Congrats it's a ZOI!
                        ZOI_calls(x,y) = true;
                        if interaction_structure.ZOI_depth(x,y)>ZOI_depth_threshold && interaction_structure.ZOI_AUC(x,y)>ZOI_AUC_upper_threshold
                            ZOI_criteria(x,y) = 3;
                        else
                            if interaction_structure.ZOI_depth(x,y)>ZOI_depth_threshold
                                ZOI_criteria = 1;
                            end
                            if interaction_structure.ZOI_AUC(x,y)>ZOI_AUC_upper_threshold
                                ZOI_criteria = 2;
                            end
                        end
                    end
                end
            end
        end
        if ZOI_criteria == 0
            error(['Failed to make ZOI call on sample ' num2str(x) ',' num2str(y)])
        end
    end
end

