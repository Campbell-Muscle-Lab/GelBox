function gui = update_display(gui,active_box)

gel_data = guidata(gui.Window);

% Display image data
set(gui.image_data_box_text, ...
    'String',printstruct(gel_data.imfinfo), ...
    'HorizontalAlignment','left');



% Display selection
if (isfield(gel_data,'box_handle'))
    n = numel(gel_data.box_handle);
    
    % Get selected box in control
    control_strings = get(gui.zoom_control,'String');
    selected_box = str2num(control_strings{ ...
                    get(gui.zoom_control,'Value')});
    
    for i=1:n
        p(i,1:4) = getPosition(gel_data.box_handle(i));
    end
    w = p(selected_box,3);
    h = p(selected_box,4);
    if ((w~=gel_data.old_width)|(h~=gel_data.old_height))
        for i=1:n
            p(i,3) = w;
            p(i,4) = h;
            setPosition(gel_data.box_handle(i),p(i,:));
        end
    end
    
    % Store data in case we need to save it
    for i=1:n
        gel_data.box_position(i,:) = getPosition(gel_data.box_handle(i));
    end
    
    % Store data for display
    d=[];
    for i=1:n
        d.box(i).fitting_mode = gel_data.fitting_mode ;
        % Extract position
        d.box(i).position = getPosition(gel_data.box_handle(i));
       
        % Label it
        set(gel_data.box_label(i),'String',sprintf('%.0f',i));
        set(gel_data.box_label(i), ...
           'Position',[d.box(i).position(1) d.box(i).position(2)]);
    
        % Calculate profile
        d.box(i).inset = imcrop(gel_data.im_data, ...
                    d.box(i).position);
        
        m = imcomplement(d.box(i).inset);
        
        x = flipud(mean(m,2));
        y = 1:size(m,1);

        x_back = linspace(x(1),x(end),numel(y));

        if (strcmp(gel_data.fitting_mode,'Double'))

            [x1,x2,x_fit] = double_gauss(y,x,x_back);
       
        else
        

        d.box(i).total_area = sum(x);
        d.box(i).background_area = sum(x_back);
        d.box(i).band_area = d.box(i).total_area - ...
                                d.box(i).background_area;
        

        end
        % Store data for later
        gel_data.box_data(i) = d.box(i);

        % Display
        if (i==selected_box)
            center_image_with_preserved_aspect_ratio( ...
                d.box(i).inset, ...
                gui.zoom_inset_axes);
       
            cla(gui.zoom_profile_axes);
            plot(gui.zoom_profile_axes,x,y,'b-');
            hold(gui.zoom_profile_axes,'on');
            plot(gui.zoom_profile_axes,x_back,y,'r-');
            
            if (strcmp(gel_data.fitting_mode,'Double'))
            fill(gui.zoom_profile_axes,x1+x_back,y,'g','FaceAlpha',0.25);
            fill(gui.zoom_profile_axes,x2+x_back,y,'m','FaceAlpha',0.25);
            plot(gui.zoom_profile_axes,x_fit+x_back,y,'k-');
            end

            xlim(gui.zoom_profile_axes,[0 max(x)]);
            xlabel(gui.zoom_profile_axes,'Intensity');
            ylabel(gui.zoom_profile_axes,'Pixels');
        end
    end
    
    gel_data.old_width = gel_data.box_position(1,3);
    gel_data.old_height = gel_data.box_position(1,4);
    
       
    % Show summary
    set(gui.selection_data_box_text, ...
        'String',printstruct(d), ...
        'HorizontalAlignment','left');
end

% Save data
guidata(gui.Window,gel_data);






