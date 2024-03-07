function new_background

no_of_panels_wide = 4;
no_of_panels_high = 1;
left_pads = 0.4 * ones(1,no_of_panels_wide);
right_pads = 0.2 * ones(1,no_of_panels_wide);
% right_pads(2) = 0.8;
panel_label_x_offset(1) = 0.2;
panel_label_x_offset(2) = 0.6;

panel_labels = '';

sp = initialise_publication_quality_figure( ...
    'no_of_panels_wide', no_of_panels_wide, ...
    'no_of_panels_high', no_of_panels_high, ...
    'top_margin', 0.15, ...
    'bottom_margin', 0.2, ...
    'right_margin', 0.005, ...
    'individual_padding', 1, ...
    'left_pads', repmat(left_pads,[1 no_of_panels_high*no_of_panels_wide]), ...
    'right_pads', repmat(right_pads,[1 no_of_panels_high*no_of_panels_wide]), ...
    'axes_padding_top', 0.2, ...
    'axes_padding_bottom',0.25, ...
    'figure_handle',1, ...
    'individual_panel_labels',repmat(panel_labels,[1 no_of_panels_wide*no_of_panels_high]));

image = imread('inset.png');
im = imcomplement(image);
% im = medfilt2(im,[3 3],'symmetric');
mean_image = mean(im,2);

mean_image = flipud(mean_image);

radius = 200;
se = strel('disk', radius);
im_close = imclose(image,se);
im_close = imcomplement(im_close);
mean_imclose = mean(im_close,2);
mean_imclose = flipud(mean_imclose);



line = linspace(mean_image(1),mean_image(end),numel(mean_image));

straight = mean_image(end)*ones(1,numel(mean_image));
subplot(sp(1))
colormap('gray')
center_image_with_preserved_aspect_ratio(image,sp(1),[])

a = mean_image;
b = line;
inter = a-b;

[xout_diff,yout_diff] = intersections(1:numel(diff(mean_image)),diff(mean_image),1:numel(diff(mean_image)),zeros(1,numel(diff(mean_image))),1);
xout_diff;
p1 = ceil(xout_diff(1)) - 1;

ix = ceil(xout_diff)-1;
d = [xout_diff];
k = 3;
[idx,C] = kmeans(d,k);

[~,c_ix] = max(C);

p2_cluster = d(idx == c_ix);
p2 = ceil(p2_cluster(1)) - 1;

background = [(mean_image(1:p1-1))'...
    linspace(mean_image(p1),mean_image(p2),(p2-p1+1))...
    (mean_image(p2+1:end))'];


subplot(sp(2))
% figure(6)
% clf
% subplot(1,3,1)
plot(mean_image,1:numel(mean_image))
ylim([1 numel(mean_image)])
hold on
plot(line,1:numel(line))
[xout,yout] = intersections(1:numel(a),a,1:numel(b),b,1);
plot(yout,xout,'ko','markersize',10)
plot(mean_image(p1),p1,'rd','markersize',10)
plot(mean_image(p2),p2,'rs','markersize',10)
plot(mean_imclose,(1:numel(mean_imclose)),'m-.')

plot(background,1:numel(background),'g-.')

xlim([0 3e4])
ylabel('Pixels')
xlabel('Density (A.U.)')
box on

% subplot(1,3,2)

subplot(sp(3))
hold on
plot(diff(mean_image),1:numel(diff(mean_image)))
plot(zeros(1,numel(diff(mean_image))),1:numel(diff(mean_image)))
plot(yout_diff,xout_diff,'ko','markersize',10)
ylim([1 numel(mean_image)])
ylabel('Pixels')
xlabel('Dens. Deriv.')
box on

subplot(sp(4))
hold on
col = parula(k);
for i = 1:k
    if i == c_ix
        c = 'k';
        m = 's';
    else
        c = col(i,:);
        m = '*';
    end
plot(d(idx==i,1),d(idx==i,1),m,'MarkerSize',12,'color',c)
end
plot(p2_cluster(1),p2_cluster(1),'s','MarkerSize',12,'color','r')
ylabel('Pixels')
xlabel('Pixels')


box on

figure_export('output_file_string','figure_new_background', ...
    'output_type','png')

figure(6)
clf
plot(mean_image,1:numel(mean_image))
ylim([1 numel(mean_image)])
hold on
plot(line,1:numel(line))
[xout,yout] = intersections(1:numel(a),a,1:numel(b),b,1);
plot(yout,xout,'ko','markersize',10)
plot(mean_image(p1),p1,'rd','markersize',10)
plot(mean_image(p2),p2,'rs','markersize',10)
plot(mean_imclose,(1:numel(mean_imclose)),'m-.')
plot(background,1:numel(background),'g-.')

legend('Raw Density','Linear Background','','','','Rolling Ball','Polyline')
ylabel('Pixels')
xlabel('Density (A.U.)')
xlim([0 3e4])

box on
figure_export('output_file_string','figure_new_background_large', ...
    'output_type','png')

end