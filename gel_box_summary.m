function gel_box_summary(summary)

[m,n] = size(summary);

total_figures = n*2;

im_handles = 1:2:total_figures;
fit_handles = 2:2:total_figures;

if mod(total_figures,6) == 0
    no_of_panels_wide = 6;
    no_of_panels_high = total_figures/no_of_panels_wide;
else
    no_of_panels_wide = 6;
    no_of_panels_high = 5;
end

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


for i = 1 : length(im_handles)
    subplot(sp(im_handles(i)))
    center_image_with_preserved_aspect_ratio( ...
        summary(i).inset,sp(im_handles(i)));

    subplot(sp(fit_handles(i)))
    plot(summary.x,summary.y,'b-');
    hold on
    plot(summary.x_fit,summary.y,'ro');

end













end