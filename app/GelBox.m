
classdef GelBox < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        AnalysisControlsPanel        matlab.ui.container.Panel
        BoxSelectionDropDown         matlab.ui.control.DropDown
        BoxSelectionDropDownLabel    matlab.ui.control.Label
        NumberofBandsDropDown        matlab.ui.control.DropDown
        NumberofBandsDropDownLabel   matlab.ui.control.Label
        GridLayout                   matlab.ui.container.GridLayout
        GelImagePanel                matlab.ui.container.Panel
        gel_image_axes               matlab.ui.control.UIAxes
        GridLayout3                  matlab.ui.container.GridLayout
        GridLayout4                  matlab.ui.container.GridLayout
        FileControlsPanel            matlab.ui.container.Panel
        OutputButton                 matlab.ui.control.Button
        SaveAnalysisButton           matlab.ui.control.Button
        LoadAnalysisButton           matlab.ui.control.Button
        InvertImageButton            matlab.ui.control.Button
        LoadImageButton              matlab.ui.control.Button
        GridLayout2                  matlab.ui.container.GridLayout
        SelectedBoxInformationPanel  matlab.ui.container.Panel
        ImageFileInformationPanel    matlab.ui.container.Panel
        BoxFitPanel                  matlab.ui.container.Panel
        box_inset                    matlab.ui.control.UIAxes
        UIAxes_2                     matlab.ui.control.UIAxes
        BoxControlsPanel             matlab.ui.container.Panel
        DeleteBoxButton              matlab.ui.control.Button
        NewBoxButton                 matlab.ui.control.Button
    end

    
    properties (Access = public)
        gel_data % Description
    end
    
    methods (Access = public)
        
        function UpdateDisplay(app)


            
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)

            addpath(genpath('../code/utilities'));
            writelines(evalc('type(mfilename(''fullpath'')+".mlapp")'),mfilename('fullpath')+".m");
%             colormap(app.gel_image_axes,'gray')
            colormap(app.UIFigure, 'gray');

        end

        % Button pushed function: LoadImageButton
        function LoadImageButtonPushed(app, event)
            [file_string,path_string]=uigetfile2( ...
                {'*.png','PNG';'*.tif','TIF'}, ...
                'Select image file');
            app.SaveAnalysisButton.Enable = 1;
            app.OutputButton.Enable = 1;
            if (path_string~=0)

                app.gel_data = [];

                app.gel_data.invert_status = 0;
                app.gel_data.image_file_string = fullfile(path_string,file_string);
                app.gel_data.im_data = imread(app.gel_data.image_file_string);
                if (ndims(app.gel_data.im_data)==3)
                    app.gel_data.im_data = rgb2gray(app.gel_data.im_data);
                end
                
                center_image_with_preserved_aspect_ratio( ...
                app.gel_data.im_data, ...
                app.gel_image_axes);

                app.gel_data.imfinfo = imfinfo(app.gel_data.image_file_string);
                app.ImInfoArea.Value = printstruct(app.gel_data.imfinfo);

            end

        end

        % Button pushed function: InvertImageButton
        function InvertImageButtonPushed(app, event)
            app.gel_data.im_data = imcomplement(app.gel_data.im_data);
            center_image_with_preserved_aspect_ratio( ...
                app.gel_data.im_data, ...
                app.gel_image_axes);
            app.gel_data.invert_status = 1;
        end

        % Button pushed function: LoadAnalysisButton
        function LoadAnalysisButtonPushed(app, event)
            % Delete any old boxes
            if (isfield(app.gel_data,'box_handle'))
                n = numel(app.gel_data.box_handle);
                for i=1:n
                    delete(app.gel_data.box_handle(i));
                end
            end
            [file_string,path_string] = uigetfile2( ...
            {'*.gdf','Gel data file'},'Select file to load prior analysis');
            
            if (path_string~=0)

            temp = load(fullfile(path_string,file_string),'-mat','save_data');
            save_data = temp.save_data;

            % Restore
            app.gel_data = [];
            app.gel_data.image_file_string = save_data.image_file_string;
            app.gel_data.im_data = save_data.im_data;
            app.gel_data.imfinfo = save_data.imfinfo;
            
            center_image_with_preserved_aspect_ratio( ...
                app.gel_data.im_data, ...
                app.gel_image_axes);
            end
            error('utku')
        end

        % Button pushed function: NewBoxButton
        function NewBoxButtonPushed(app, event)
            if (~isfield(app.gel_data,'box_handle'))
                app.DeleteBoxButton.Enable = 1;
                n=1;
                app.gel_data.box_handle(n) = drawrectangle(app.gel_image_axes);
                p = app.gel_data.box_handle(n).Position;
                app.gel_data.old_width = p(3);
                app.gel_data.old_height = p(4);
            else
                n = 1 + numel(app.gel_data.box_handle);
                p = app.gel_data.box_handle(n-1).Position;

                app.gel_data.box_handle(n) = images.roi.Rectangle(app.gel_image_axes, ...
                    'Position',p + [20,0,0,0]);
                for i=1:(n-1)
                    app.gel_data.box_handle(i).InteractionsAllowed = 'none';
                end
            end

            addlistener(app.gel_data.box_handle(n),"MovingROI", ...
                @(src,evt) new_box_position(evt));

            % Set color to last box
            app.gel_data.box_handle(n).Color = [0 1 0];
            app.gel_data.box_handle(n).FaceAlpha = 0;
            for i=1:(n-1)
                app.gel_data.box_handle(i).Color = [1 0 0];
            end

            % Add in a label
            p = app.gel_data.box_handle(n).Position;
            app.gel_data.box_label(n) = text(app.gel_image_axes, ...
                p(1)+p(3),p(2)-50,sprintf('%.0f',n));

            % Update zoom control
            for i=1:n
                control_strings{i}=sprintf('%.0f',i);
            end

            app.BoxSelectionDropDown.Items = control_strings;
            app.BoxSelectionDropDown.Value = control_strings{n};


            function new_box_position(evt);
                if (isfield(app.gel_data,'box_position'))
                    box_position = app.gel_data.box_position;
                    [r,c]=size(box_position);
                    if (r>=n)&(~isequal(box_position(n,:),evt.CurrentPosition))
%                         update_display(gui,n);
                    end
                else
%                     update_display(gui,n);
                end
            end

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            colormap(app.UIFigure, 'parula');
            app.UIFigure.Position = [100 100 1119 608];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.WindowState = 'maximized';

            % Create BoxControlsPanel
            app.BoxControlsPanel = uipanel(app.UIFigure);
            app.BoxControlsPanel.Title = 'Box Controls';
            app.BoxControlsPanel.Position = [191 -140 234 64];

            % Create NewBoxButton
            app.NewBoxButton = uibutton(app.BoxControlsPanel, 'push');
            app.NewBoxButton.ButtonPushedFcn = createCallbackFcn(app, @NewBoxButtonPushed, true);
            app.NewBoxButton.Position = [9 13 100 22];
            app.NewBoxButton.Text = 'New Box';

            % Create DeleteBoxButton
            app.DeleteBoxButton = uibutton(app.BoxControlsPanel, 'push');
            app.DeleteBoxButton.Enable = 'off';
            app.DeleteBoxButton.Position = [123 13 100 22];
            app.DeleteBoxButton.Text = 'Delete Box';

            % Create BoxFitPanel
            app.BoxFitPanel = uipanel(app.UIFigure);
            app.BoxFitPanel.Title = 'Box Fit (?)';
            app.BoxFitPanel.Position = [526 -658 572 381];

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.BoxFitPanel);
            app.UIAxes_2.Position = [165 -1 394 330];

            % Create box_inset
            app.box_inset = uiaxes(app.BoxFitPanel);
            app.box_inset.Position = [14 11 134 330];

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.RowHeight = {'1x', '2x'};

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.GridLayout);
            app.GridLayout2.RowHeight = {'1x'};
            app.GridLayout2.Layout.Row = 1;
            app.GridLayout2.Layout.Column = 2;

            % Create ImageFileInformationPanel
            app.ImageFileInformationPanel = uipanel(app.GridLayout2);
            app.ImageFileInformationPanel.Title = 'Image File Information';
            app.ImageFileInformationPanel.Layout.Row = 1;
            app.ImageFileInformationPanel.Layout.Column = 1;

            % Create SelectedBoxInformationPanel
            app.SelectedBoxInformationPanel = uipanel(app.GridLayout2);
            app.SelectedBoxInformationPanel.Title = 'Selected Box Information';
            app.SelectedBoxInformationPanel.Layout.Row = 1;
            app.SelectedBoxInformationPanel.Layout.Column = 2;

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.GridLayout);
            app.GridLayout3.ColumnWidth = {'1x', '3x'};
            app.GridLayout3.RowHeight = {'1x'};
            app.GridLayout3.Layout.Row = 1;
            app.GridLayout3.Layout.Column = 1;

            % Create FileControlsPanel
            app.FileControlsPanel = uipanel(app.GridLayout3);
            app.FileControlsPanel.Title = 'File Controls';
            app.FileControlsPanel.Layout.Row = 1;
            app.FileControlsPanel.Layout.Column = 1;

            % Create LoadImageButton
            app.LoadImageButton = uibutton(app.FileControlsPanel, 'push');
            app.LoadImageButton.ButtonPushedFcn = createCallbackFcn(app, @LoadImageButtonPushed, true);
            app.LoadImageButton.Position = [19 128 90 22];
            app.LoadImageButton.Text = 'Load Image';

            % Create InvertImageButton
            app.InvertImageButton = uibutton(app.FileControlsPanel, 'push');
            app.InvertImageButton.ButtonPushedFcn = createCallbackFcn(app, @InvertImageButtonPushed, true);
            app.InvertImageButton.Position = [19 96 90 22];
            app.InvertImageButton.Text = 'Invert Image';

            % Create LoadAnalysisButton
            app.LoadAnalysisButton = uibutton(app.FileControlsPanel, 'push');
            app.LoadAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @LoadAnalysisButtonPushed, true);
            app.LoadAnalysisButton.Position = [19 64 90 22];
            app.LoadAnalysisButton.Text = 'Load Analysis';

            % Create SaveAnalysisButton
            app.SaveAnalysisButton = uibutton(app.FileControlsPanel, 'push');
            app.SaveAnalysisButton.Enable = 'off';
            app.SaveAnalysisButton.Position = [19 33 90 23];
            app.SaveAnalysisButton.Text = 'Save Analysis';

            % Create OutputButton
            app.OutputButton = uibutton(app.FileControlsPanel, 'push');
            app.OutputButton.Enable = 'off';
            app.OutputButton.Position = [19 2 90 22];
            app.OutputButton.Text = 'Output';

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.GridLayout3);
            app.GridLayout4.RowHeight = {'1x', '1x', '1x'};
            app.GridLayout4.Layout.Row = 1;
            app.GridLayout4.Layout.Column = 2;

            % Create GelImagePanel
            app.GelImagePanel = uipanel(app.GridLayout);
            app.GelImagePanel.Title = 'Gel Image';
            app.GelImagePanel.Layout.Row = 2;
            app.GelImagePanel.Layout.Column = 1;

            % Create gel_image_axes
            app.gel_image_axes = uiaxes(app.GelImagePanel);
            app.gel_image_axes.Visible = 'off';
            app.gel_image_axes.Position = [23 24 468 330];

            % Create AnalysisControlsPanel
            app.AnalysisControlsPanel = uipanel(app.UIFigure);
            app.AnalysisControlsPanel.Title = 'Analysis Controls';
            app.AnalysisControlsPanel.Position = [-236 -312 386 173];

            % Create NumberofBandsDropDownLabel
            app.NumberofBandsDropDownLabel = uilabel(app.AnalysisControlsPanel);
            app.NumberofBandsDropDownLabel.Position = [9 117 99 22];
            app.NumberofBandsDropDownLabel.Text = 'Number of Bands';

            % Create NumberofBandsDropDown
            app.NumberofBandsDropDown = uidropdown(app.AnalysisControlsPanel);
            app.NumberofBandsDropDown.Items = {'1', '2', '3'};
            app.NumberofBandsDropDown.Position = [123 117 100 22];
            app.NumberofBandsDropDown.Value = '1';

            % Create BoxSelectionDropDownLabel
            app.BoxSelectionDropDownLabel = uilabel(app.AnalysisControlsPanel);
            app.BoxSelectionDropDownLabel.Position = [9 81 98 22];
            app.BoxSelectionDropDownLabel.Text = 'Box Selection';

            % Create BoxSelectionDropDown
            app.BoxSelectionDropDown = uidropdown(app.AnalysisControlsPanel);
            app.BoxSelectionDropDown.Items = {};
            app.BoxSelectionDropDown.Placeholder = 'No data';
            app.BoxSelectionDropDown.Position = [123 81 100 22];
            app.BoxSelectionDropDown.Value = {};

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GelBox

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end

