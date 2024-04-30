function [center,radius,metric,intense_med,intense_mean,intense_std,intense_serr] = imfind_radialintensity(img,object_range,scan_range)
    
[H,W] = size(img);
max_datapoints = 2*ceil(2*pi*scan_range(2)); %safety factor of 2 in case between pixels
intensities = nan(max_datapoints,scan_range(2));

[centers, radii, metrics] = imfindcircles(img,object_range,'Sensitivity',0.99);
if isempty(metrics)
    metric=0;
    center=[0 0];
    radius=0;
    intense_mean=nan(1,scan_range(2));
    intense_med=nan(1,scan_range(2));
    intense_std=nan(1,scan_range(2));
    intense_serr=nan(1,scan_range(2));
    return
else
    metric = metrics(1);
    center = centers(1,:);
    radius = radii(1);
end
intense_n = zeros(1,scan_range(2));
% intense_S1 = zeros(1,scan_range(2));
% intense_S2 = zeros(1,scan_range(2));
for x=1:W
    for y=1:H
        hdist = int16(hypot(x-centers(1,1),y - centers(1,2)));
        if (hdist <= scan_range(2)) & (hdist>0)
            intense_n(hdist) = intense_n(hdist) + 1;
%             intense_S1(hdist) =  intense_S1(hdist) + single(img(y,x));            
%             intense_S2(hdist) = intense_S2(hdist) + single(img(y,x)).^2;
            intensities(intense_n(hdist),hdist) = img(y,x); %rows columns
        end
    end
end
intense_med=single(nanmedian(intensities));
intense_mean=nanmean(single(intensities));
intense_std=nanstd(single(intensities));
intense_serr=intense_std./sqrt(intense_n);

% intense_mean=single(intense_S1)./intense_n; 
% intense_std=(sqrt(single(intense_S2./intense_n)-single((intense_S1./intense_n)).^2));
