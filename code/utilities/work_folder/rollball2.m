function rollball2


inset = imread('math_image.png');
m = 4;

figure(1)
clf

colormap('gray')
a = subplot(1,m,1);
imagesc(inset)
set(a,'YDir','Reverse');

u = subplot(1,m,2);
se = strel('disk',100);
inset_d = imdilate(inset,se);
inset_d = imerode(inset_d,se);

imagesc(inset_d)
set(u,'YDir','Reverse');

u = subplot(1,m,3);
hold on
plot(65535-mean(inset,2),1:size(inset,2))
plot(65535-mean(inset_d,2),1:size(inset,2),'o')
set(u,'YDir','Reverse');

end