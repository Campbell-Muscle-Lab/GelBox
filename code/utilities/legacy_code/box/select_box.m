function gui = new_box(gui);

gel_data = guidata(gui.Window)

if (~isfield(gel_data,'box'))
    n=0;
else
    n = numel(gel_data.box);
end

n=n+1;

gel_data.box(n) = imrect(gui.gel_axes);
gel_data.box(n).position = getPosition(gel_data.box(n));
addNewPositionCallback(gel_data.box(n),@new_box_position);

guidata(gui.Window,gel_data)

    % Nested function
    function new_box_position(pos);
        gel_data = guidata(gui.Window);
        update_display(gui);
    end
        
end
        


    

