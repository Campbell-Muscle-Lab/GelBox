function gui = update_display(gui,active_box)

gel_data = guidata(gui.Window);

temp_strings = get(gui.fitting_mode,'String');
gel_data.fitting_mode = temp_strings{get(gui.fitting_mode,'Value')};

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

        num_of_bands = str2double(gel_data.fitting_mode);
        [x_bands,x_fit] = fit_gaussian(y,x,x_back,num_of_bands);
        
        d.box(i).total_area = sum(x);
        d.box(i).background_area = sum(x_back);
        % Check if band number is changed
        if size(x_bands,1) ~= num_of_bands
            num_of_bands = size(x_bands,1);
        end
        
        for j = 1 : num_of_bands
            d.box(i).band_area(j) = trapz(y,x_bands(j,:));
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
            plot(gui.zoom_profile_axes,x_fit+x_back,y,'k-');
            color = lines(num_of_bands);

            for j = 1 : num_of_bands

                fill(gui.zoom_profile_axes,x_back+x_bands(j,:),y,color(j,:),'FaceAlpha',0.25)
            
            end
            plot(gui.zoom_profile_axes,x_back,y,'r-');
            x_limit = max(max(x),max(x_fit));
            xlim(gui.zoom_profile_axes,[0 x_limit+10]);
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






