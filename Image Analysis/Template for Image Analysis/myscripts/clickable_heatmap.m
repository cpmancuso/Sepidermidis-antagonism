function clickable_heatmap(imgprefix,all_adj_centers,all_radii,norm_x,norm_int1,norm_int2,ZOI_call,thresh, headers)

figure(2)
clf(2)

im = imagesc(ZOI_call');
im.ButtonDownFcn = @mouse_click;
hold on

xlabel('Reciever Lawn')
ylabel('Producer Spot')

[num ~] = size(ZOI_call');

xrange = [1 num];
dx = diff(xrange)/(num-1);
xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,num+1);
hm = mesh(xg,xg,zeros(num+1));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';


    function mouse_click(src,eventData)
        crop_dim=load('params.mat', 'crop_dim').crop_dim;
        crop_radius=load('params.mat', 'crop_radius').crop_radius;
        int_radius=load('params.mat', 'int_radius').int_radius;

        % get coordinates of click 
        coords = eventData.IntersectionPoint;

        i=round(coords(1));
        n=round(coords(2));
        img_i = headers.imagenum(i);
        
        disp([i,n])
        %load image from masked image saved in folder
        f=figure(3);
        clf(f)
        f.Position=[0.2000 206.6000 1.4856e+03 642.4000];
        subplot(1,2,1)
        filename = ['Masked Images\',imgprefix,sprintf('%02d',img_i),'_rect_masked.jpg'];
        new_img = imcrop(imread(filename),[all_adj_centers(i,n,1)-crop_radius,all_adj_centers(i,n,2)-crop_radius,2*crop_radius,2*crop_radius]);
        imshow(new_img)
        hold on
        viscircles(gca,[crop_radius+1,crop_radius+1],all_radii(i,n),'Color','r','LineWidth',0.1)

        subplot(1,2,2)

        plot(reshape(norm_x(i,n,:),[],int_radius),reshape(norm_int1(i,n,:),[],int_radius),'b')
        hold on
        plot(reshape(norm_x(i,n,:),[],int_radius),reshape(norm_int2(i,n,:),[],int_radius),'k')
        plot([-100,100],[0,0],'r-','LineWidth',0.5)
        plot([-100,100],[-thresh,-thresh],'r--','LineWidth',0.5)
        plot([0,0],[-50,75],'r-','LineWidth',0.5)
        ylim([-50,75])
        title(['Spot ' num2str(n) ' on Lawn ' num2str(i)])
        legend('Circular Background Correction','Square Background Correction','','','Threshold')
        xlabel('Pixels from Colony Edge')
        ylabel('Normalized Intensity')
    end
end



