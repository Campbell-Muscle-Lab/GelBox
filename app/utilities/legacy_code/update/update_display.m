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
        p(i,1:4) = gel_data.box_handle(i).Position;
    end
    w = p(selected_box,3);
    h = p(selected_box,4);
    if ((w~=gel_data.old_width)|(h~=gel_data.old_height))
        for i=1:n
            p(i,3) = w;
            p(i,4) = h;
            gel_data.box_handle(i).Position = p(i,:);
        end
    end

    % Store data in case we need to save it
    for i=1:n
        gel_data.box_position(i,:) = gel_data.box_handle(i).Position;
    end

    % Store data for display
    d=[];
    for i=1:n
        d.box(i).fitting_mode = gel_data.fitting_mode ;
        % Extract position
        d.box(i).position = gel_data.box_handle(i).Position;

        % Label it
        set(gel_data.box_label(i),'String',sprintf('%.0f',i));
        set(gel_data.box_label(i), ...
            'Position',[d.box(i).position(1)+d.box(i).position(3) ...
            d.box(i).position(2)-50]);

        % Calculate profile
        d.box(i).inset = imcrop(gel_data.im_data, ...
            d.box(i).position);

        m = imcomplement(d.box(i).inset);

        x = flipud(mean(m,2));
        y = 1:size(m,1);

        x_back = linspace(x(1),x(end),numel(y));

        num_of_bands = str2double(gel_data.fitting_mode);
        if gel_data.loaded_analysis || gel_data.parameters_updated || gel_data.par_update(i)
        elseif gel_data.new_box || gel_data.mode_updated || gel_data.par_est_na
            [par_est,par_con] = estimate_fitting_parameters(y,x, ...
                x_back,num_of_bands);
            fnames = fieldnames(par_est);
            for k = 1:numel(fnames)
                gel_data.par_est(i).(fnames{k}) = [];
                gel_data.par_est(i).(fnames{k}) = ...
                    par_est.(fnames{k});
                gel_data.par_con(i).(fnames{k}) = [];
                gel_data.par_con(i).(fnames{k}) = ...
                    par_con.(fnames{k});
            end
        end
        switch num_of_bands
            case 2
                [x_bands,x_fit,r_squared] = ...
                fit_2gaussian(y,x,x_back,i);
            case 3
                [x_bands,x_fit,r_squared] = ...
                fit_3gaussian(gel_data,y,x,x_back,i);
        end
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
        gel_data.summary(i).x = x;
        gel_data.summary(i).y = y;
        gel_data.summary(i).x_fit = x_fit;
        gel_data.summary(i).x_back = x_back;

        if num_of_bands == 2
            [~,peak_band_1] = max(x_bands(1,:));    
            [~,peak_band_2] = max(x_bands(2,:));

            if peak_band_2 > peak_band_1

                gel_data.summary(i).bottom = d.box(i).band_area(1);
                gel_data.summary(i).top = d.box(i).band_area(2);
                gel_data.summary(i).band_1 = x_bands(1,:);
                gel_data.summary(i).band_2 = x_bands(2,:);
            else
                gel_data.summary(i).bottom = d.box(i).band_area(2);
                gel_data.summary(i).top = d.box(i).band_area(1);
                gel_data.summary(i).band_1 = x_bands(2,:);
                gel_data.summary(i).band_2 = x_bands(1,:);
            end
        elseif num_of_bands == 3
            [~,peak_band_1] = max(x_bands(1,:));    
            [~,peak_band_2] = max(x_bands(2,:));
            [~,peak_band_3] = max(x_bands(3,:));
            
            peak_band = [peak_band_1 peak_band_2 peak_band_3];

            %sort
            [~,sort_ix] = sort(peak_band);
            gel_data.summary(i).bottom = d.box(i).band_area(sort_ix(1));
            gel_data.summary(i).middle = d.box(i).band_area(sort_ix(2));
            gel_data.summary(i).top = d.box(i).band_area(sort_ix(3));
            gel_data.summary(i).band_1 = x_bands(sort_ix(1),:);
            gel_data.summary(i).band_2 = x_bands(sort_ix(2),:);
            gel_data.summary(i).band_3 = x_bands(sort_ix(3),:);
        else
            gel_data.summary(i).band_1 = x_bands(1,:);
        end


        gel_data.summary(i).inset = d.box(i).inset;
        gel_data.summary(i).r_squared = r_squared;


        % Display
        if (i==selected_box)
            center_image_with_preserved_aspect_ratio( ...
                d.box(i).inset, ...
                gui.zoom_inset_axes);

            cla(gui.zoom_profile_axes);
            cla(gui.zoom_profile_axes_corrected);
            plot(gui.zoom_profile_axes,x,y,'k-');
            hold(gui.zoom_profile_axes,'on');
            plot(gui.zoom_profile_axes,x_fit+x_back,y,':r','LineWidth',2);
            
            switch num_of_bands
                case 2
                    color = {'r','b'};
                case 3
                    color = {'r','b','g'};
            end
            for j = 1 : num_of_bands

                fill(gui.zoom_profile_axes,x_back+x_bands(j,:),y,color{j},'FaceAlpha',0.25)

            end
            plot(gui.zoom_profile_axes,x_back,y,'m-');
            x_limit = max(max(x),max(x_fit+x_back));
            xlabel(gui.zoom_profile_axes,'Optical Density');
            ylabel(gui.zoom_profile_axes,'Pixels');
            x_t_end = ceil(max(x)/50)*50;
            x_t_mid = round(x_t_end/2);
            try
                x_ticks = [0 x_t_mid x_t_end];
                xticks(gui.zoom_profile_axes,x_ticks)
                xlim(gui.zoom_profile_axes,[0 x_t_end])
            end
            ylim(gui.zoom_profile_axes,[1 y(end)])
            plot(gui.zoom_profile_axes_corrected,x-x_back',y,'k', ...
                'LineWidth',1.5)
            hold(gui.zoom_profile_axes_corrected,'on')

            t = sprintf('Box %i r^2 = %.3f',i,r_squared);
            title(gui.zoom_profile_axes_corrected,t)
            area_1 = trapz(y,x_bands(1,:));
            area_2 = trapz(y,x_bands(2,:));
            plot(gui.zoom_profile_axes_corrected,x_bands(1,:),y,'ro');
            plot(gui.zoom_profile_axes_corrected,x_bands(2,:),y,'bo');
            try
            area_3 = trapz(y,x_bands(3,:));
            plot(gui.zoom_profile_axes_corrected,x_bands(3,:),y,'go');
            end
            switch num_of_bands
                case 2
                    area_tot = area_1 + area_2;
                    str1 = sprintf('Red: %.2f\n Blue: %.2f', ...
                        area_1/area_tot, area_2/area_tot);
                case 3
                    area_tot = area_1 + area_2 + area_3;
                    str1 = sprintf('Red: %.2f\n Blue: %.2f\n Green: %.2f', ...
                        area_1/area_tot, area_2/area_tot, area_3/area_tot);
            end

            plot(gui.zoom_profile_axes_corrected,x_fit,y,':r', ...
                'LineWidth',1.5)

            xlabel(gui.zoom_profile_axes_corrected,'Back. Corr. Density')
            ylabel(gui.zoom_profile_axes_corrected,'Pixels')
            ylim(gui.zoom_profile_axes_corrected,[1 y(end)])
            xL=xlim(gui.zoom_profile_axes_corrected);
            yL=ylim(gui.zoom_profile_axes_corrected);
            text(gui.zoom_profile_axes_corrected,0.99*xL(2),0.99*yL(2),str1, ...
                'HorizontalAlignment','right','VerticalAlignment','top')

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






