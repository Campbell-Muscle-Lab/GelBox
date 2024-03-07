function rollball


inset = imread('math_image.png');

im = imcomplement(inset);
mean_inset = mean(im,2);
% mean_inset = flipud(mean_inset);

mean_inset_ui = uint16(mean_inset);
radius = 150;
se = strel('disk',radius);
im_closing = imdilate(inset,se);
corr_1D = imopen(mean_inset_ui,se);
corr = imsubtract(inset,imcomplement(im_closing));
corr_mean = mean(imcomplement(corr),2);
m = 4

figure(1)
colormap('gray')
clf
h = subplot(1,m,1);
imagesc((inset))
title('Math Image')
set(h,'YDir','Reverse');
h2 = subplot(1,m,2);
plot(mean_inset,1:numel(mean_inset))

set(h2,'YDir','Reverse');
subplot(1,m,3)
imagesc(corr)
h3 = subplot(1,m,4)
plot(mean_inset,1:numel(mean_inset))
hold on
plot(corr_mean,1:numel(mean_inset))
plot(corr_1D,1:numel(corr_1D))
set(h3,'YDir','Reverse');





return












mean_fiji = 65535-flipud(mean_fiji);
mean_fiji_p = mean_fiji;

mean_inset = mean(inset,2);
mean_inset = 65535 - flipud(mean_inset);

% se = offsetstrel('ball', 80,80,0);
% mean_im_close = imclose(mean_inset,se);

% mean_inset = 

radius = [ 80];
colors = jet(5);
% figure(2)
% clf
% hold on
% plot(1:numel(mean_inset),mean_inset,'color','k')
% plot(1:numel(mean_fiji_p),mean_fiji_p,'color','r')
% plot(1:numel(mean_fiji_p),-mean_fiji_p+mean_inset)
% 


for i = 1
se = strel('disk', radius(i));
im_close = imclose(mean(inset,2),se);

mean_imclose = mean(im_close,2);
mean_imclose = flipud(mean_imclose);

figure(2)
clf
hold on 
% plot(1:numel(mean_imclose),65535 - mean_imclose,'color',colors(i,:),'linewidth',1.5)
plot(1:numel(mean_imclose),65535 - mean_imclose,'color','g','linewidth',1.5)
plot(1:numel(mean_imclose),65535-mean_inset - mean_imclose,'o','color','r','linewidth',1.5)
xlabel('Pixels')
ylabel('Density (A.U.)')

end


xlim([1, numel(mean_inset)])








end