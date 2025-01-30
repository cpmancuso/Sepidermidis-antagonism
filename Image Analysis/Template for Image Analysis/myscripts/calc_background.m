function [background_circle,background_square] = calc_background(img,back_radius)
% Takes a masked image centered on the center of a masked spot, 
% calculates circular background and square background

[H,W] = size(img);
center = [(W+1)/2 (H+1)/2]; %assumes image is centered

img_background = single(img);
img_background(img_background==0) = NaN; %Mask gets saved as 0. Assumes darkest areas have non-zero int, but this is not ideal


[xx,yy] = ndgrid((1:H)-center(2),(1:W)-center(1));
roi = (xx.^2 + yy.^2)<back_radius^2;
background_circle=nanmedian(img_background(roi));
background_square=nanmedian(nanmedian(img_background));
