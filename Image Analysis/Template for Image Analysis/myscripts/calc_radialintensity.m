function [intense_med,intense_mean,intense_std,intense_serr] = calc_radialintensity(img,int_radius)
% Takes an image centered on the center of a spot, calculates radial intensity
[H,W] = size(img);
max_datapoints = 2*ceil(2*pi*int_radius); %safety factor of 2 in case between pixels
intensities = nan(max_datapoints,int_radius);
center = [(W+1)/2 (H+1)/2]; %assumes image is centered
intense_n = zeros(1,int_radius);

for w=1:W
    for h=1:H
        hdist = int16(hypot(w-center(1),h - center(2)));
        if (hdist <= int_radius) && (hdist>0)
            intense_n(hdist) = intense_n(hdist) + 1; %number of pixels at this radius
            if img(h,w) ~= 0 %Mask gets saved as 0. Assumes darkest areas have non-zero int, but this is not ideal
                intensities(intense_n(hdist),hdist) = img(h,w); %rows columns
            end
        end
    end
end
intense_med=single(nanmedian(intensities));
intense_mean=nanmean(single(intensities));
intense_std=nanstd(single(intensities));
intense_serr=intense_std./sqrt(intense_n);
