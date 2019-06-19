function gui = zoom_control_update(~,~,gui)

gel_data = guidata(gui.Window);

% Get selected box in control
control_strings = get(gui.zoom_control,'String');
selected_box = str2num(control_strings{ ...
                    get(gui.zoom_control,'Value')});
n = numel(control_strings);                

for i=1:n
    if (i~=selected_box)
        setColor(gel_data.box_handle(i),[1 0 0]);
        setResizable(gel_data.box_handle(i),false);
    else
        setColor(gel_data.box_handle(i),[0 1 0]);
        setResizable(gel_data.box_handle(i),true);
    end
end

update_display(gui,selected_box);

