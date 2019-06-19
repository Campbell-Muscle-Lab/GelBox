function gui = new_box(gui);

gel_data = guidata(gui.Window)

if (~isfield(gel_data,'box_handle'))
    n=1;
    gel_data.box_handle(n) = imrect(gui.gel_axes);
    p = getPosition(gel_data.box_handle(n));
    gel_data.old_width = p(3);
    gel_data.old_height = p(4);
else
    n = 1 + numel(gel_data.box_handle);
    p = getPosition(gel_data.box_handle(n-1));
    gel_data.box_handle(n) = imrect(gui.gel_axes,p);
    setPosition(gel_data.box_handle(n),p+[20 0 0 0]);
    for i=1:(n-1)
        setResizable(gel_data.box_handle(i),false);
    end
end

guidata(gui.Window,gel_data)
addNewPositionCallback(gel_data.box_handle(n),@new_box_position);

% Set color to last box
setColor(gel_data.box_handle(n),[0 1 0]);
for i=1:(n-1)
    setColor(gel_data.box_handle(i),[1 0 0]);
end

% Add in a label
p = getPosition(gel_data.box_handle(n));
gel_data.box_label(n) = text(p(1),p(2),sprintf('%.0f',n));

% Update zoom control
for i=1:n
    control_strings{i}=sprintf('%.0f',i);
end
set(gui.zoom_control,'String',control_strings);
set(gui.zoom_control,'Value',n);

guidata(gui.Window,gel_data)

update_display(gui,n);

    % Nested function
    function new_box_position(pos);
        gel_data = guidata(gui.Window);
        if (isfield(gel_data,'box_position'))
            box_position = gel_data.box_position;
            [r,c]=size(box_position);
            if (r>=n)&(~isequal(box_position(n,:),pos))
                update_display(gui,n);
            end
        else
            update_display(gui,n);
        end
    end
        
end
        


    

