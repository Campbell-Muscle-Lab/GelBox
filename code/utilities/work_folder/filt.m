function filt



im = imread('inset2.png');


im_f = medfilt2(im,[3 3],'symmetric');

im_f2 = medfilt2(im,[10 10],'symmetric');

figure(1)
colormap('gray')
clf
subplot(1,3,1)
imshow(im)
title('Original')
subplot(1,3,2)
imshow(im_f)
title('Filter size of 3')

subplot(1,3,3)
imshow(im_f2)
title('Filter size of 10')

% figure(2)
% clf
% imshowpair(im,im_f,"montage")



figure_export('output_file_string','medfilt', ...
    'output_type','png')




end