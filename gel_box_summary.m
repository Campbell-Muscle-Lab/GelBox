function gel_box_summary(summary,path_string,file_string)

[m,n] = size(summary);

path_string = regexprep(path_string,'excel_files','output_figures')
output_file_string = fullfile(path_string,file_string);
output_file_string = erase(output_file_string,'.xlsx');
fig_title = erase(file_string,'.xlsx');
fig_title = regexprep(fig_title,'_',' ');
output_file_string_b = sprintf('%s_boxes',output_file_string);
output_file_string_r = sprintf('%s_r_squared',output_file_string);

total_figures = n*2;

im_handles = 1:2:total_figures;
fit_handles = 2:2:total_figures;
no_of_panels_wide = 4;
no_of_panels_high = (total_figures + ...
    mod(total_figures,no_of_panels_wide))/no_of_panels_wide;


left_pads = 0.35 * ones(1,no_of_panels_wide);
right_pads = 0.15 * ones(1,no_of_panels_wide);

sp = initialise_publication_quality_figure( ...
    'no_of_panels_wide', no_of_panels_wide, ...
    'no_of_panels_high', no_of_panels_high, ...
    'top_margin', 0.2, ...
    'bottom_margin', 0, ...
    'right_margin', 1, ...
    'individual_padding', 1, ...
    'left_pads', repmat(left_pads,[1 no_of_panels_high]), ...
    'right_pads', repmat(right_pads,[1 no_of_panels_high]), ...
    'axes_padding_top', 0.2, ...
    'axes_padding_bottom',0.5, ...
    'panel_label_font_size', 0, ...
    'figure_handle',19);

numel(sp)
for i = 1 : length(im_handles)
    h = subplot(sp(im_handles(i)));
    colormap(h,"gray")
    t = sprintf('%s Box %i',fig_title,i);
    title(t)
    center_image_with_preserved_aspect_ratio( ...
        summary(i).inset,sp(im_handles(i)));

   
    subplot(sp(fit_handles(i)))
    plot(summary(i).x,summary(i).y,'k-');
    hold on
    fill(sp(fit_handles(i)),summary(i).x_back+summary(i).band_1, ...
        summary(i).y,'b','FaceAlpha',0.1)
    fill(sp(fit_handles(i)),summary(i).x_back+summary(i).band_2, ...
        summary(i).y,'y','FaceAlpha',0.1)
    plot(summary(i).x_back,summary(i).y, ...
        'mo','MarkerSize',1);
    plot(summary(i).x_fit+summary(i).x_back,summary(i).y, ...
        'ro','MarkerSize',1);
     xlabel(sp(fit_handles(i)),'Intensity');
     ylabel(sp(fit_handles(i)),'Pixels');
    t = sprintf('r^2 = %.3f',summary(i).r_squared);
    title(t)
    str1 = sprintf('Top: %.3f\n Bottom: %.3f', summary(i).top, summary(i).bottom);

    xL=xlim;
    yL=ylim;
    text(0.99*xL(2),0.99*yL(2),str1,'HorizontalAlignment','right','VerticalAlignment','top')


end


% improve_axes

figure_export('output_file_string', output_file_string_b,'output_type', 'png');

figure(20)
cla
for i = 1 : length(im_handles)
    plot(i, summary(i).r_squared,'bo')
    hold on
    xlabel('Box No')
    ylabel('r squared')

end

xlim([0.9 n*1.1])
xticks([1:n])
ylim([0 1.05])
figure_export('output_file_string', output_file_string_r,'output_type', 'png');



end