function corners = calibrate_platemap(calibration_figure)

title('Carefully click the center of the 4 corner wells')

xvals = [];
yvals = [];

for num = 1:4
    [x,y] = ginput(1);
    xvals(num) = x;
    yvals(num) = y;
    hold on
    plot(x,y,'o','MarkerSize', 4, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
end

corners = [xvals; yvals];
close(calibration_figure)
