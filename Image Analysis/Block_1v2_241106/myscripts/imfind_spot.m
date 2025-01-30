function [center,radius,metric] = imfind_spot(crop_img,crop_dim,object_range)  

[centers, radii, metrics] = imfindcircles(crop_img,object_range,'Sensitivity',1);
if isempty(metrics) % no spot detected
    metric=0;
    center=[round(crop_dim/2) round(crop_dim/2)];
    radius=0;
    intense_mean=nan(1,object_range(2));
    intense_med=nan(1,object_range(2));
    intense_std=nan(1,object_range(2));
    intense_serr=nan(1,object_range(2));
    return
else % use best spot guess
    metric = metrics(1);
    center = centers(1,:);
    radius = radii(1);
end

