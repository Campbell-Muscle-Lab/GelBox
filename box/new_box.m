function gui = new_box(gui);

gel_data = guidata(gui.Window)

if (~isfield(gel_data,'box_handle'))
    n=1;
    gel_data.box_handle(n) = drawrectangle(gui.gel_axes);
    p = gel_data.box_handle(n).Position;
    gel_data.old_width = p(3);
    gel_data.old_height = p(4);
else
    n = 1 + numel(gel_data.box_handle);
    p = gel_data.box_handle(n-1).Position;
    
    gel_data.box_handle(n) = images.roi.Rectangle(gui.gel_axes,'Position',p + [20,0,0,0]);
    for i=1:(n-1)
        gel_data.box_handle(i).InteractionsAllowed = 'none';
    end
end

guidata(gui.Window,gel_data)
addlistener(gel_data.box_handle(n),"ROIMoved",@(src,evt) new_box_position(evt));

% Set color to last box
gel_data.box_handle(n).Color = [0 1 0];
gel_data.box_handle(n).FaceAlpha = 0;
for i=1:(n-1)
    gel_data.box_handle(i).Color = [1 0 0];
end

% Add in a label
p = gel_data.box_handle(n).Position;
gel_data.box_label(n) = text(p(1)+p(3),p(2)-50,sprintf('%.0f',n));

% Update zoom control
for i=1:n
    control_strings{i}=sprintf('%.0f',i);
end
set(gui.zoom_control,'String',control_strings);
set(gui.zoom_control,'Value',n);

guidata(gui.Window,gel_data)

update_display(gui,n);

% Nested function
    function new_box_position(evt);
        gel_data = guidata(gui.Window);
        if (isfield(gel_data,'box_position'))
            box_position = gel_data.box_position;
            [r,c]=size(box_position);
            if (r>=n)&(~isequal(box_position(n,:),evt.CurrentPosition))
                update_display(gui,n);
            end
        else
            update_display(gui,n);
        end
    end

end





