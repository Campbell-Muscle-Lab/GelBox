function gui = create_interface(version_string)
% Create the interface

gui.version_string = version_string;

% Open a window
screen_size = get(0,'ScreenSize');
gui.Window = figure( ...
    'Name',gui.version_string, ...
    'NumberTitle','off', ...
    'Toolbar','figure', ...
    'Menubar','none', ...
    'Resize','off', ...
    'HandleVisibility','on', ...
    'Position',[0.1 0.1 0.8 0.8].*screen_size([3 4 3 4]));
set(gui.Window,'Colormap',gray);

% Arrange the main window
gui.main_layout = uix.HBox( ...
    'Parent',gui.Window, ...
    'Spacing',3);

% Divide into panels
gui.gel_panel = uix.Panel( ...
    'Parent',gui.main_layout);
gui.zoom_panel = uix.Panel( ...
    'Parent',gui.main_layout);
gui.data_panel = uix.Panel( ...
    'Parent',gui.main_layout);
set(gui.main_layout,'Widths',[-3 -1 -2]);

% Add in gel axes
gui.gel_box = uix.HBox( ...
    'Parent',gui.gel_panel);
gui.gel_axes = axes( ...
    'Parent',uicontainer('Parent',gui.gel_box), ...
    'Visible','off');

% Add in zoom display
gui.zoom_box = uix.VBox( ...
    'Parent',gui.zoom_panel);
gui.zoom_control_panel = uix.Panel( ...
    'Parent',gui.zoom_box);
gui.zoom_inset_panel = uix.Panel( ...
    'Parent',gui.zoom_box);
gui.zoom_profile_panel = uix.Panel( ...
    'Parent',gui.zoom_box);
set(gui.zoom_box,'Heights',[-1 -3 -3]);

% Add in zoom features
gui.zoom_control_box = uix.HBox( ...
    'Parent',gui.zoom_control_panel);
gui.zoom_control_text = uicontrol( ...
    'parent',gui.zoom_control_box, ...
    'style','text', ...
    'string','Adjustable box');
gui.zoom_control = uicontrol( ...
    'parent',gui.zoom_control_box, ...
    'style','popupmenu', ...
    'string',{'No data'});
gui.zoom_inset_axes = axes( ...
    'Parent',gui.zoom_inset_panel, ...
    'Visible','off');
gui.zoom_profile_axes = axes( ...
    'Parent',gui.zoom_profile_panel, ...
    'Visible','off');

% Add in data display panels
gui.data_box = uix.VBox( ...
    'Parent',gui.data_panel);
gui.image_info_panel = uix.Panel( ...
    'Parent',gui.data_box);
gui.selection_data_panel = uix.Panel( ...
    'Parent',gui.data_box);
set(gui.data_box,'Heights',[-1 -3]);
gui.image_data_box_text = uicontrol( ...
    'style','edit', ...
    'parent',gui.image_info_panel, ...
    'enable','inactive', ...
    'min',0, ...
    'max',2);
gui.selection_data_box_text = uicontrol( ...
    'style','edit', ...
    'parent',gui.selection_data_panel, ...
    'enable','inactive', ...
    'min',0, ...
    'max',2);


% Add in menus
gui.file_menu = uimenu(gui.Window,'Label','File');
uimenu(gui.file_menu, ...
    'Label','Load image', ...
    'Callback',{@call_load_image,gui});
uimenu(gui.file_menu, ...
    'Label','Save analysis', ...
    'Separator','on', ...
    'Callback',{@call_save_analysis,gui});
uimenu(gui.file_menu, ...
    'Label','Load analysis', ...
    'Callback',{@call_load_analysis,gui});
uimenu(gui.file_menu, ...
    'Label','Output results to Excel', ...
    'Separator','on', ...
    'Callback',{@call_output_results,gui});


gui.edit_menu = uimenu(gui.Window,'Label','Edit');
uimenu(gui.edit_menu, ...
    'Label','New box', ...
    'Callback',{@call_new_box,gui}, ...
    'Accelerator','n');

gui.tools_menu = uimenu(gui.Window,'Label','Tools');
uimenu(gui.tools_menu, ...
    'Label','Invert image', ...
    'Callback',{@call_invert_image,gui});


% Add in more callbacks - these have to come late to include entries
set(gui.zoom_control,'callback',{@zoom_control_update,gui});

% Nested functions
    function gui = call_load_image(~,~,gui)
        
        [file_string,path_string]=uigetfile2( ...
            {'*.png','PNG';'*.tif','TIF'}, ...
            'Select image file');
        
        if (path_string~=0)
            
            gel_data = [];
            
            gel_data.invert_status = 0;
            
            gel_data.image_file_string = fullfile(path_string,file_string);
            gel_data.im_data = imread(gel_data.image_file_string);
            if (ndims(gel_data.im_data)==3)
                gel_data.im_data = rgb2gray(gel_data.im_data);
            end

            center_image_with_preserved_aspect_ratio( ...
                gel_data.im_data, ...
                gui.gel_axes);
            
            gel_data.imfinfo = imfinfo(gel_data.image_file_string);
            
            guidata(gui.gel_axes,gel_data);
            
            update_display(gui)
        end
    end

    function gui = call_new_box(~,~,gui)
        new_box(gui);
    end

    function gui = call_invert_image(~,~,gui)
        gel_data = guidata(gui.Window);
        gel_data.im_data = imcomplement(gel_data.im_data);
        center_image_with_preserved_aspect_ratio( ...
                gel_data.im_data, ...
                gui.gel_axes);
        gel_data.invert_status = 1;
        guidata(gui.Window,gel_data);
    end
        
    function gui = call_save_analysis(~,~,gui)
        gel_data = guidata(gui.Window);
        
        save_data.image_file_string = gel_data.image_file_string;
        save_data.box_position = gel_data.box_position;     
        save_data.im_data = gel_data.im_data;
        save_data.imfinfo = gel_data.imfinfo;
        
        [file_string,path_string] = uiputfile2( ...
            {'*.gdf','Gel data file'},'Select file to save analysis');

        if (path_string~=0)
            save(fullfile(path_string,file_string),'save_data');
           
            msgbox(sprintf('Current analysis saved to %s',file_string), ...
                'Analysis saved');
        end

    end

    function gui = call_load_analysis(~,~,gui)

        gel_data = guidata(gui.Window);
        
        % Delete any old boxes
        if (isfield(gel_data,'box_handle'))
            n = numel(gel_data.box_handle);
            for i=1:n
                delete(gel_data.box_handle(i));
            end
        end
        
        [file_string,path_string] = uigetfile2( ...
            {'*.gdf','Gel data file'},'Select file to load prior analysis');
        
        if (path_string~=0)


            temp = load(fullfile(path_string,file_string),'-mat','save_data');
            save_data = temp.save_data;

            % Restore
            gel_data = [];
            gel_data.image_file_string = save_data.image_file_string;
            gel_data.im_data = save_data.im_data;
            gel_data.imfinfo = save_data.imfinfo;

            center_image_with_preserved_aspect_ratio( ...
                gel_data.im_data, ...
                gui.gel_axes);

            n=size(save_data.box_position,1);
            control_strings = [];
            for i=1:n
                gel_data.box_handle(i) = ...
                    imrect(gui.gel_axes,save_data.box_position(i,:));
                control_strings{i} = sprintf('%.0f',i);
            end
            set(gui.zoom_control,'String',control_strings);
            set(gui.zoom_control,'Value',1);

            gel_data.box_label=[];

            for i=1:n
                if (i~=1)
                    setColor(gel_data.box_handle(i),[1 0 0]);
                    setResizable(gel_data.box_handle(i),false);
                else
                    setColor(gel_data.box_handle(i),[0 1 0]);
                    setResizable(gel_data.box_handle(i),true);
                end

                p = getPosition(gel_data.box_handle(i))
                gel_data.box_label(i) = text(p(1),p(2),sprintf('%.0f',i), ...
                    'Parent',gui.gel_axes);

                gel_data.old_width = p(3);
                gel_data.old_height = p(4);


                i=i
                addNewPositionCallback(gel_data.box_handle(i),@new_box_position2);
            end

            % Need this to make labels
            drawnow;
        end

        guidata(gui.gel_axes,gel_data);
        update_display(gui)
        
            % Nested function
            function new_box_position2(pos);
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

    function gui = call_output_results(~,~,gui)
        
        gel_data = guidata(gui.Window)
        
        % Save data as structure for output
        d = [];
        d.image_file{1} = gel_data.image_file_string;
        n = numel(gel_data.box_handle);
        for i=1:n
            d.band(i) = i;
            d.total_area(i) = gel_data.box_data(i).total_area;
            d.background_area(i) = gel_data.box_data(i).background_area;
            d.band_area(i) = gel_data.box_data(i).band_area;
            d.band_left(i) = gel_data.box_data(i).position(1);
            d.band_top(i) = gel_data.box_data(i).position(2);
            d.band_width(i) = gel_data.box_data(i).position(3);
            d.band_height(i) = gel_data.box_data(i).position(4);
        end
        
       [file_string,path_string] = uiputfile2( ...
            {'*.xlsx','Excel file'},'Select file for results');
        
        if (path_string~=0)
            write_structure_to_excel( ...
                'filename',fullfile(path_string,file_string), ...
                'structure',d);
            
            msgbox(sprintf('Data written to %s',file_string), ...
                'Data saved');
        end
    end
end
        








