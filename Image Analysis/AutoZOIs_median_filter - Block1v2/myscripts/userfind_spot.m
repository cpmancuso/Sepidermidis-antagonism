function [center,radius,metric] = userfind_spot(img,object_range)

figure(3)
imshow(img)
roi = drawcircle;

metric = 1;
center = roi.Center;
radius = roi.Radius;
