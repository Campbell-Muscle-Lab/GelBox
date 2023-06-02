function center_image_with_preserved_aspect_ratio(im,axis_handle)

if (nargin==1)
    axis_handle = gca;
end

old_axis_units=get(axis_handle,'Units');
axes_pos=get(axis_handle,'Position');
set(axis_handle,'Units','Pixels');
axes_pos=get(axis_handle,'Position');
axes_x=axes_pos(3);
axes_y=axes_pos(4);

im_x=size(im,2);
im_y=size(im,1);

axes_aspect_ratio=axes_x/axes_y;
im_aspect_ratio=im_x/im_y;

% Display the image
imagesc(im,'Parent',axis_handle);

if (axes_aspect_ratio>im_aspect_ratio)
    % We have to pad the image horizontally so that it has the same
    % aspect ratio as the axes
    
    new_im_x=im_y*axes_aspect_ratio;
    pad_x=(new_im_x-im_x)/2;

    xlim(axis_handle,[-pad_x im_x+pad_x]);
    ylim(axis_handle,[1 im_y]);
    
else
    % We have to pad the image vertically so that it has the same
    % aspect ratio as the axes

    new_im_y=im_x/axes_aspect_ratio;
    pad_y=(new_im_y-im_y)/2;
    
    xlim(axis_handle,[1 im_x]);
    ylim(axis_handle,[-pad_y im_y+pad_y]);
end

set(axis_handle,'YDir','Reverse');
set(axis_handle,'Visible','off');
set(findall(axis_handle, 'type', 'text'), 'visible', 'on')
set(axis_handle,'Units',old_axis_units);