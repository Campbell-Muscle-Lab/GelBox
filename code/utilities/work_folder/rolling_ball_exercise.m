function rolling_ball(image, radius)


im = imread('inset.png');
im_rolled = imread('inset_rolled_80.png');
% im = uint8(im / 256);
mean_int = mean(imcomplement(im),2);
mean_int = flipud(mean_int);
mean_int_r = mean(imcomplement(im_rolled),2);
mean_int_r = flipud(mean_int_r);
benchmark = mean_int-mean_int_r;
figure(1)
clf
plot(mean_int)
hold on
plot(mean_int_r)
% x = load('ecg.mat');
% 
% ecg_10 = load('ecg_10.mat');
% ecg_10 = ecg_10.ecg_10;
% 
% ecg_80 = load('ecg_80.mat');
% ecg_80 = ecg_80.ecg_80;

im_c = imcomplement(im);

rad = 80;
se = offsetstrel('ball',rad,rad);
% se = strel('disk',rad);
% figure(88)
% imshow(se.Neighborhood);
% return
% utku = imtophat(im_2,se);
% imwrite(utku,'example_rolled_50.png')
% figure(4)
% imshow(imcomplement(utku))
eroded_im_c = imclose(im_c,se);
% x = 65535 - back;
% eroded_im_c = imadd(im_c,x);

% h = (2^16-1) - eroded_im_c;
% eroded_im_c = imadd(eroded_im_c,h);
mean_int_err = mean((eroded_im_c),2);
mean_int_err = flipud(mean_int_err);

im_corrected = im_c - eroded_im_c;

figure(2)
clf
% subplot(2,6,1)
% imshow(im_c)
% title('Original')
% 
% subplot(2,6,2)
% imshow(imcomplement(im_c));
% title({'"Complement"', 'Original'})
% 
% subplot(2,6,3)
% imshow(eroded_im_c)
% title('Erroded')
% 
% subplot(2,6,4)
% imshow(imcomplement(eroded_im_c))
% title({'"Complement"', 'Erroded'})
% 
% subplot(2,6,5)
% imshow(im_corrected)
% title('Corrected')
% 
% subplot(2,6,6)
% imshow(imcomplement(im_corrected));
% title({'"Complement"', 'Corrected'})
% 
% % montage({im_c,im_c-eroded_im_c,imcomplement(eroded_im_c)},'size',[1,3])
% 
mean_int_cor = mean((im_c-eroded_im_c),2);
mean_int_cor = flipud(mean_int_cor);
% subplot(2,1,2)
plot(mean_int)
hold on
plot(mean_int_cor)
plot(mean_int_err)
plot(mean_int-mean_int_r,'o')
plot(mean_int_r,'s')

xlim([1 numel(mean_int)])
% xlabel('Pixels')
% ylabel({'Optical','Density (A.U.)'})
% legend('Mean: Complement Original','Mean: Complement Corrected','Original - Erroded')
return
figure_export('output_file_string','20 radius', ...
    'output_type','png')
return
prof = readmatrix('profile.txt');
prof = [zeros(1,10)'; prof; zeros(1,10)'];

% ecg = x.utku;
lin_back = linspace(prof(1),prof(end),numel(prof))';


radius = 20; 
se = strel('ball', radius, radius, 0);

% filtered_ecg_erode = imerode(ecg, se) + radius;
% filtered_ecg_tophat = imtophat(ecg,se);
filtered_prof_erode = imerode(prof, se);
filtered_prof_tophat = imtophat(prof, se);

double_prof_tophat = imtophat(prof-lin_back, se);


figure(1)
clf
hold on

plot(prof,'linewidth',1.75)
plot(filtered_prof_erode,'linewidth',1.75)
% plot(prof-filtered_prof_tophat,'linewidth',1.75)
% plot(prof-lin_back,'linewidth',1.75)
% plot(double_prof_tophat,'linewidth',1.75)


% legend('Profile','Rolling: Imerode','Rolling: Imtophat','Linear Correction','Lin+Roll')

% % plot(ecg_10)
% plot(ecg_80,'o','color','k')
% plot(ecg-filtered_ecg_erode,'linewidth',1.5)
% plot(filtered_ecg_tophat,'linewidth',1.5)

% plot(prof)
% plot(prof-filtered_prof)
% legend('original','me','scikit10','scikit80')
% ylabel('Voltage')
% xlabel('Time')
% xlim([1 numel(ecg)])
box on 
end
