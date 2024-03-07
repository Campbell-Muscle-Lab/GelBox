function math_image

addpath(genpath('../../MATLAB_Utilities'))
no_of_panels_wide = 4;
no_of_panels_high = 3;
left_pads = 0.4 * ones(1,no_of_panels_wide);
right_pads = 0.2 * ones(1,no_of_panels_wide);
right_pads(2) = 0.8;
panel_label_x_offset(1) = 0.2;
panel_label_x_offset(2) = 0.6;

panel_labels = '';

sp = initialise_publication_quality_figure( ...
    'no_of_panels_wide', no_of_panels_wide, ...
    'no_of_panels_high', no_of_panels_high, ...
    'top_margin', 0.15, ...
    'bottom_margin', 0, ...
    'right_margin', 1, ...
    'individual_padding', 1, ...
    'left_pads', repmat(left_pads,[1 no_of_panels_high*no_of_panels_wide]), ...
    'right_pads', repmat(right_pads,[1 no_of_panels_high*no_of_panels_wide]), ...
    'axes_padding_top', 0.2, ...
    'axes_padding_bottom',0.25, ...
    'figure_handle',1, ...
    'individual_panel_labels',repmat(panel_labels,[1 no_of_panels_wide*no_of_panels_high]));


utku = load('matlab.mat');
bands = utku.x_bands;
trapz(bands(:,1));
trapz(bands(:,2));
n = 150;
x = linspace(1, 300,n);

xx = linspace(1,150,n);

x = x';

y1 = (10*x)*1;
y2 = 0.1*x.^2;
y3 = 0.0005*x.^3;

x0 = [100 200];
gamma = [0.001 0.001];
skew1 = [0 0];
A = [8000 10000];

g1  = skewed_Gaussian(x,x0(1),gamma(1),A(1),skew1(1));
g2  = skewed_Gaussian(x,x0(2),gamma(2),A(2),skew1(2));

p1 = y1 + g1 + g2;
p2 = y2 + g1 + g2;
p3 = y3 + g1 + g2;

subplot(sp(1))
plot(y1,xx)
title('Line')
subplot(sp(2))
hold on
f1 = plot(g1,xx);
f2 = plot(g2,xx);
str1 = sprintf('Area = %.4f',simps(xx,g1)/(simps(xx,g1)+simps(xx,g2)));
str2 = sprintf('Area = %.4f',simps(xx,g2)/(simps(xx,g1)+simps(xx,g2)));

figs = [f1,f2];
leg_labels = {str1, str2};

legendflex(figs, leg_labels, ...
    'xscale',0.45, ...
    'anchor',{'ne','ne'}, ...
    'buffer',[70 -5], ...
    'padding',[1 1 2], ...
    'FontSize',7, ...
    'text_y_padding', -1);
title('Gaussians')
subplot(sp(3))
plot(p1,xx)
title('Line + Gaussians')

math_image1 = repmat(p1,[1 n]);
math_image1 = uint16(math_image1);

subplot(sp(4))
colormap('gray')
imagesc(imcomplement(math_image1))
title('Pseudo Image')

imwrite(imcomplement(math_image1),'math_image_1st.png')



% 


subplot(sp(5))
plot(y2,xx)
title('x^{2} Based')
subplot(sp(6))
hold on
f1 = plot(g1,xx);
f2 = plot(g2,xx);
str1 = sprintf('Area = %.4f',simps(xx,g1)/(simps(xx,g1)+simps(xx,g2)));
str2 = sprintf('Area = %.4f',simps(xx,g2)/(simps(xx,g1)+simps(xx,g2)));

figs = [f1,f2];
leg_labels = {str1, str2};

legendflex(figs, leg_labels, ...
    'xscale',0.45, ...
    'anchor',{'ne','ne'}, ...
    'buffer',[70 -5], ...
    'padding',[1 1 2], ...
    'FontSize',7, ...
    'text_y_padding', -1);
title('Gaussians')
subplot(sp(7))
plot(p2,xx)
title('x^{2} Based + Gaussians')

math_image2 = repmat(p2,[1 n]);
math_image2 = uint16(math_image2);

subplot(sp(8))
colormap('gray')
imagesc(imcomplement(math_image2))
set(sp(1:12),'YDir','Reverse');
title('Pseudo Image')
imwrite(imcomplement(math_image2),'math_image_2nd.png')


subplot(sp(9))
plot(y3,xx)
title('x^{3} Based')
subplot(sp(10))
hold on
f1 = plot(g1,xx);
f2 = plot(g2,xx);
str1 = sprintf('Area = %.4f',simps(xx,g1)/(simps(xx,g1)+simps(xx,g2)));
str2 = sprintf('Area = %.4f',simps(xx,g2)/(simps(xx,g1)+simps(xx,g2)));

figs = [f1,f2];
leg_labels = {str1, str2};

legendflex(figs, leg_labels, ...
    'xscale',0.45, ...
    'anchor',{'ne','ne'}, ...
    'buffer',[70 -5], ...
    'padding',[1 1 2], ...
    'FontSize',7, ...
    'text_y_padding', -1);
title('Gaussians')
subplot(sp(11))
plot(p3,xx)
title('x^{3} Based + Gaussians')

math_image3 = repmat(p3,[1 n]);
math_image3 = uint16(math_image3);

subplot(sp(12))
colormap('gray')
imagesc(imcomplement(math_image3))
set(sp(1:12),'YDir','Reverse');
title('Pseudo Image')
imwrite(imcomplement(math_image3),'math_image_3rd.png')

pad = zeros(150,25);

pseudo_gel = [pad math_image1 pad math_image2 pad math_image3 pad];
imwrite(imcomplement(pseudo_gel),'pseudo_gel.png')

figure_export('output_file_string','figure_math_image', ...
    'output_type','png')

    function y=skewed_Gaussian(x,x0,gamma,A,skew1)
        offset = zeros(length(x),1);
        offset((x-x0)>0) = skew1*(x((x-x0)>0)-x0);
        y=  A*exp(-gamma*(((x-x0)+offset).^2));
    end




end