function [center,radius,metric] = imfind_spot(img,object_range)  

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

