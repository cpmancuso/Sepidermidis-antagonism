clc 
clear all
close all

for f=1:2
    % Open an image
    % [file, imagePath] = uigetfile('C:\Users\cmanc\Dropbox\Lieberman Lab\Personal lab notebooks\Chris Mancuso\Data\2023_06_05\*.jpg');  % Replace with the actual path to your image
    [file, imagePath] = uigetfile('C:\Users\cmanc\Dropbox\Lieberman Lab\Personal lab notebooks\Chris Mancuso\Data\2023_06_28\*.jpg');  % Replace with the actual path to your image
    img = imread([imagePath,file]);
    
    % Display the image
    imshow(img);
    title('Click on the image to select a crop area');
    
    % Wait for user to click on the image
    [x, y] = ginput(1);  % Get the coordinates of the clicked point
    
    % Close the image window
    close(gcf);
    
    % Define the size of the crop
    cropSize = 400;
    
    % Calculate the crop region
    xmin = max(1, round(x - cropSize/2));
    ymin = max(1, round(y - cropSize/2));
    width = min(size(img, 2) - xmin, cropSize);
    height = min(size(img, 1) - ymin, cropSize);
    
    % Crop the image
    croppedImage = img(ymin:ymin+height-1, xmin:xmin+width-1, :);
    
    imshow(croppedImage)
    % Save the cropped image
    imwrite(croppedImage, ['Selected Images\257vs227_' num2str(f) '.png']);  % Replace with the desired save path
end