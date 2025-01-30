function [center,radius,metric] = userfind_spot(img,crop_dim,object_range)

figure(3)
clf(3)
imshow(img)
hold on
[rows,cols,~] = size(img);

quiver(0.7*rows, 0.9*cols, rows/5, -cols/5,'r-', 'LineWidth', 1,'MaxHeadSize',10);
title('Click and drag to exactly surround spot')
roi = drawcircle;

metric = 1;
center = roi.Center;
radius = roi.Radius;
