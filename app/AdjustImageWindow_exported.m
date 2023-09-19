classdef AdjustImageWindow_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        AdjustImageUIFigure            matlab.ui.Figure
        BrightnessandContrastPanel     matlab.ui.container.Panel
        ContrastLamp                   matlab.ui.control.Lamp
        MaximizeContrastButton         matlab.ui.control.Button
        ContrastUpperBoundSlider       matlab.ui.control.Slider
        ContrastUpperBoundLabel        matlab.ui.control.Label
        ContrastLowerBoundSlider       matlab.ui.control.Slider
        ContrastLowerBoundSliderLabel  matlab.ui.control.Label
        AdjustBrightnessSlider         matlab.ui.control.Slider
        AdjustBrightnessSliderLabel    matlab.ui.control.Label
        adjusted_image_hist            matlab.ui.control.UIAxes
        original_image_hist            matlab.ui.control.UIAxes
        AdjustedImagePanel             matlab.ui.container.Panel
        RevertChangesButton            matlab.ui.control.Button
        LoadAdjustedImageButton        matlab.ui.control.Button
        adjusted_image_axis            matlab.ui.control.UIAxes
        OriginalImagePanel             matlab.ui.container.Panel
        CropImageButton                matlab.ui.control.Button
        original_image_axis            matlab.ui.control.UIAxes
    end


    properties (Access = private)
        GelBoxApp % Description
        n
    end

    properties (Access = public)
        original_image
        im_crop_box
        adjusted_image = []
        adjusted_image_2 = []
        normalized_original_image % Description
        cropped_image % Description
        crop_pos % Description
        bit_d_val % Description
        max_contrast = false
    end

    methods (Access = public)

        function UpdateAdjustedImageDisplay(app,image)

            cla(app.adjusted_image_axis);
            cla(app.adjusted_image_hist);

            center_image_with_preserved_aspect_ratio( ...
                image, ...
                app.adjusted_image_axis,[0 1]);

            [counts, edges] = histcounts(image,100);
            bar(app.adjusted_image_hist,edges(1:(end-1)), counts);
            

        end
        
        function LoadImageAdjustments(app)
            
            app.original_image = app.GelBoxApp.gel_data.original_image;
            
            bit_d = class(app.original_image);

            switch bit_d
                case 'uint8'
                    app.bit_d_val = 2^8;
                    app.original_image = double(app.original_image)./(app.bit_d_val);
                case 'uint16'
                    app.bit_d_val = 2^16;
                    app.original_image = double(app.original_image)./(app.bit_d_val);
            end


            [counts, edges] = histcounts(app.original_image,100);
            bar(app.original_image_hist,edges(1:(end-1)), counts);
            center_image_with_preserved_aspect_ratio( ...
                app.original_image, ...
                app.original_image_axis,[]);
            
            app.adjusted_image = app.GelBoxApp.gel_data.adjusted_image;
            bit_d = class(app.adjusted_image);

            switch bit_d
                case 'uint8'
                    app.bit_d_val = 2^8;
                    app.adjusted_image = double(app.adjusted_image)./(app.bit_d_val);
                case 'uint16'
                    app.bit_d_val = 2^16;
                    app.adjusted_image = double(app.adjusted_image)./(app.bit_d_val);
            end

            UpdateAdjustedImageDisplay(app,app.adjusted_image)
            
            app.AdjustBrightnessSlider.Value = app.GelBoxApp.gel_data.image_adjustments.brightness;
            app.ContrastLowerBoundSlider.Value = app.GelBoxApp.gel_data.image_adjustments.contrast_lower;
            app.ContrastUpperBoundSlider.Value = app.GelBoxApp.gel_data.image_adjustments.contrast_upper;
            
            if app.GelBoxApp.gel_data.image_adjustments.max_contrast
                app.ContrastLamp.Enable = 'on';
            end

        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, caller)
            movegui(app.AdjustImageUIFigure,'center')
            colormap(app.AdjustImageUIFigure, 'gray');
            xlim(app.adjusted_image_hist,[0 1]);
            xlim(app.adjusted_image_hist,[0 1]);
            cla(app.original_image_hist);
            app.GelBoxApp = caller;
            
            if isfield(app.GelBoxApp.gel_data,'image_adjustments')
                LoadImageAdjustments(app)
            else
                app.original_image = caller.gel_data.im_data;

                center_image_with_preserved_aspect_ratio( ...
                    app.original_image, ...
                    app.original_image_axis,[]);

                bit_d = class(app.original_image);

                switch bit_d
                    case 'uint8'
                        app.bit_d_val = 2^8;
                        app.original_image = double(app.original_image)./(app.bit_d_val);
                    case 'uint16'
                        app.bit_d_val = 2^16;
                        app.original_image = double(app.original_image)./(app.bit_d_val);
                end


                [counts, edges] = histcounts(app.original_image,100);
                bar(app.original_image_hist,edges(1:(end-1)), counts);
                UpdateAdjustedImageDisplay(app,app.original_image)
            end

        end

        % Button pushed function: CropImageButton
        function CropImageButtonPushed(app, event)
            app.im_crop_box = drawrectangle(app.original_image_axis);
            app.im_crop_box.FaceAlpha = 0;
            app.crop_pos = app.im_crop_box.Position;

            app.cropped_image = imcrop(app.original_image,app.crop_pos);

            addlistener(app.im_crop_box,"MovingROI", ...
                @(src,evt) crop_box_position(evt));

            UpdateAdjustedImageDisplay(app,app.cropped_image)

            function crop_box_position(evt);
                app.cropped_image = [];
                app.crop_pos = app.im_crop_box.Position;
                app.cropped_image = imcrop(app.original_image,app.crop_pos);
                UpdateAdjustedImageDisplay(app,app.cropped_image);
            end

        end

        % Value changing function: AdjustBrightnessSlider
        function AdjustBrightnessSliderValueChanging(app, event)
            changingValue = event.Value;
            if app.ContrastLowerBoundSlider.Value == 0 && app.ContrastUpperBoundSlider.Value == 1
                if changingValue < 0
                    app.adjusted_image = imadjust(app.cropped_image,[0 1], [0 1+changingValue]);
                    UpdateAdjustedImageDisplay(app,app.adjusted_image);
                else
                    app.adjusted_image = imadjust(app.cropped_image,[0 1], [changingValue 1]);
                    UpdateAdjustedImageDisplay(app,app.adjusted_image);
                end
            else
                if changingValue < 0
                    app.adjusted_image_2 = imadjust(app.adjusted_image,[0 1], [0 1+changingValue]);
                    UpdateAdjustedImageDisplay(app,app.adjusted_image_2);
                else
                    app.adjusted_image_2 = imadjust(app.adjusted_image,[0 1], [changingValue 1]);
                    UpdateAdjustedImageDisplay(app,app.adjusted_image_2);
                end
            end
        end

        % Button pushed function: MaximizeContrastButton
        function MaximizeContrastButtonPushed(app, event)
            if isempty(app.adjusted_image)
                app.adjusted_image = imadjust(app.cropped_image);
                UpdateAdjustedImageDisplay(app,app.adjusted_image);
            else
                app.adjusted_image_2 = imadjust(app.adjusted_image);
                UpdateAdjustedImageDisplay(app,app.adjusted_image_2);
            end
            
            app.max_contrast = true;
            app.ContrastLamp.Enable = 'on';
            
        end

        % Value changing function: ContrastLowerBoundSlider
        function ContrastLowerBoundSliderValueChanging(app, event)
            changingValue = event.Value;
            if app.AdjustBrightnessSlider.Value == 0
               app.adjusted_image = imadjust(app.cropped_image,[changingValue app.ContrastUpperBoundSlider.Value],[]);
               UpdateAdjustedImageDisplay(app,app.adjusted_image);
            else
               app.adjusted_image_2 = imadjust(app.adjusted_image,[changingValue app.ContrastUpperBoundSlider.Value],[]);
               UpdateAdjustedImageDisplay(app,app.adjusted_image_2);
            end
        end

        % Value changing function: ContrastUpperBoundSlider
        function ContrastUpperBoundSliderValueChanging(app, event)
            changingValue = event.Value;
            if app.AdjustBrightnessSlider.Value == 0
                app.adjusted_image = imadjust(app.cropped_image,[app.ContrastLowerBoundSlider.Value changingValue],[]);
                UpdateAdjustedImageDisplay(app,app.adjusted_image);
            else
                app.adjusted_image_2 = imadjust(app.adjusted_image,[app.ContrastLowerBoundSlider.Value changingValue],[]);
                UpdateAdjustedImageDisplay(app,app.adjusted_image_2);
            end
        end

        % Button pushed function: RevertChangesButton
        function RevertChangesButtonPushed(app, event)
            app.adjusted_image = [];
            app.adjusted_image_2 = [];
            app.ContrastLowerBoundSlider.Value = 0;
            app.ContrastUpperBoundSlider.Value = 1;
            app.AdjustBrightnessSlider.Value = 0;
            UpdateAdjustedImageDisplay(app,app.cropped_image)
        end

        % Button pushed function: LoadAdjustedImageButton
        function LoadAdjustedImageButtonPushed(app, event)
           
            app.GelBoxApp.gel_data.image_adjustments.crop_pos = app.crop_pos;
            app.GelBoxApp.gel_data.image_adjustments.brightness = app.AdjustBrightnessSlider.Value;
            app.GelBoxApp.gel_data.image_adjustments.contrast_lower = app.ContrastLowerBoundSlider.Value;
            app.GelBoxApp.gel_data.image_adjustments.contrast_upper = app.ContrastUpperBoundSlider.Value;
            app.GelBoxApp.gel_data.image_adjustments.max_contrast = app.max_contrast;
     
            if isempty(app.cropped_image)
                app.GelBoxApp.gel_data.adjusted_image = [];
            elseif isempty(app.adjusted_image)
                app.GelBoxApp.gel_data.adjusted_image = app.cropped_image.*app.bit_d_val;
            elseif isempty(app.adjusted_image_2)
                app.GelBoxApp.gel_data.adjusted_image = app.adjusted_image.*app.bit_d_val;
            else
                app.GelBoxApp.gel_data.adjusted_image = app.adjusted_image_2.*app.bit_d_val;
            end
            
            switch app.bit_d_val
                case 2^8
                    app.GelBoxApp.gel_data.adjusted_image = uint8(app.GelBoxApp.gel_data.adjusted_image);
                case 2^16
                    app.GelBoxApp.gel_data.adjusted_image = uint16(app.GelBoxApp.gel_data.adjusted_image);
            end
            ImageAdjusted(app.GelBoxApp)
            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create AdjustImageUIFigure and hide until all components are created
            app.AdjustImageUIFigure = uifigure('Visible', 'off');
            app.AdjustImageUIFigure.Position = [100 100 1722 585];
            app.AdjustImageUIFigure.Name = 'Adjust Image';

            % Create OriginalImagePanel
            app.OriginalImagePanel = uipanel(app.AdjustImageUIFigure);
            app.OriginalImagePanel.Title = 'Original Image';
            app.OriginalImagePanel.Position = [13 17 664 557];

            % Create original_image_axis
            app.original_image_axis = uiaxes(app.OriginalImagePanel);
            app.original_image_axis.XTick = [];
            app.original_image_axis.YTick = [];
            app.original_image_axis.Box = 'on';
            app.original_image_axis.Position = [6 11 637 475];

            % Create CropImageButton
            app.CropImageButton = uibutton(app.OriginalImagePanel, 'push');
            app.CropImageButton.ButtonPushedFcn = createCallbackFcn(app, @CropImageButtonPushed, true);
            app.CropImageButton.Position = [15 501 100 22];
            app.CropImageButton.Text = 'Crop Image';

            % Create AdjustedImagePanel
            app.AdjustedImagePanel = uipanel(app.AdjustImageUIFigure);
            app.AdjustedImagePanel.Title = 'Adjusted Image';
            app.AdjustedImagePanel.Position = [1045 17 664 557];

            % Create adjusted_image_axis
            app.adjusted_image_axis = uiaxes(app.AdjustedImagePanel);
            app.adjusted_image_axis.XTick = [];
            app.adjusted_image_axis.YTick = [];
            app.adjusted_image_axis.Box = 'on';
            app.adjusted_image_axis.Position = [6 11 637 475];

            % Create LoadAdjustedImageButton
            app.LoadAdjustedImageButton = uibutton(app.AdjustedImagePanel, 'push');
            app.LoadAdjustedImageButton.ButtonPushedFcn = createCallbackFcn(app, @LoadAdjustedImageButtonPushed, true);
            app.LoadAdjustedImageButton.Position = [23 497 128 22];
            app.LoadAdjustedImageButton.Text = 'Load Adjusted Image';

            % Create RevertChangesButton
            app.RevertChangesButton = uibutton(app.AdjustedImagePanel, 'push');
            app.RevertChangesButton.ButtonPushedFcn = createCallbackFcn(app, @RevertChangesButtonPushed, true);
            app.RevertChangesButton.Position = [167 497 102 22];
            app.RevertChangesButton.Text = 'Revert Changes';

            % Create BrightnessandContrastPanel
            app.BrightnessandContrastPanel = uipanel(app.AdjustImageUIFigure);
            app.BrightnessandContrastPanel.Title = 'Brightness and Contrast';
            app.BrightnessandContrastPanel.Position = [688 17 348 557];

            % Create original_image_hist
            app.original_image_hist = uiaxes(app.BrightnessandContrastPanel);
            title(app.original_image_hist, 'Original Image Histogram')
            xlabel(app.original_image_hist, 'Pixel Intensity')
            ylabel(app.original_image_hist, 'Number of Pixels')
            app.original_image_hist.Box = 'on';
            app.original_image_hist.Position = [24 370 300 162];

            % Create adjusted_image_hist
            app.adjusted_image_hist = uiaxes(app.BrightnessandContrastPanel);
            title(app.adjusted_image_hist, 'Adjusted Image Histogram')
            xlabel(app.adjusted_image_hist, 'Pixel Intensity')
            ylabel(app.adjusted_image_hist, 'Number of Pixels')
            app.adjusted_image_hist.Box = 'on';
            app.adjusted_image_hist.Position = [24 197 300 162];

            % Create AdjustBrightnessSliderLabel
            app.AdjustBrightnessSliderLabel = uilabel(app.BrightnessandContrastPanel);
            app.AdjustBrightnessSliderLabel.HorizontalAlignment = 'center';
            app.AdjustBrightnessSliderLabel.WordWrap = 'on';
            app.AdjustBrightnessSliderLabel.Position = [41 153 57 30];
            app.AdjustBrightnessSliderLabel.Text = 'Adjust Brightness';

            % Create AdjustBrightnessSlider
            app.AdjustBrightnessSlider = uislider(app.BrightnessandContrastPanel);
            app.AdjustBrightnessSlider.Limits = [-1 1];
            app.AdjustBrightnessSlider.MajorTicks = [-1 -0.6 -0.2 0 0.2 0.6 1];
            app.AdjustBrightnessSlider.MajorTickLabels = {'Darken', '', '', '', '', '', 'Brighten'};
            app.AdjustBrightnessSlider.ValueChangingFcn = createCallbackFcn(app, @AdjustBrightnessSliderValueChanging, true);
            app.AdjustBrightnessSlider.Position = [130 174 150 3];

            % Create ContrastLowerBoundSliderLabel
            app.ContrastLowerBoundSliderLabel = uilabel(app.BrightnessandContrastPanel);
            app.ContrastLowerBoundSliderLabel.HorizontalAlignment = 'center';
            app.ContrastLowerBoundSliderLabel.WordWrap = 'on';
            app.ContrastLowerBoundSliderLabel.Position = [33 111 72 29];
            app.ContrastLowerBoundSliderLabel.Text = 'Contrast Lower Bound';

            % Create ContrastLowerBoundSlider
            app.ContrastLowerBoundSlider = uislider(app.BrightnessandContrastPanel);
            app.ContrastLowerBoundSlider.Limits = [0 1];
            app.ContrastLowerBoundSlider.ValueChangingFcn = createCallbackFcn(app, @ContrastLowerBoundSliderValueChanging, true);
            app.ContrastLowerBoundSlider.Position = [130 127 150 3];

            % Create ContrastUpperBoundLabel
            app.ContrastUpperBoundLabel = uilabel(app.BrightnessandContrastPanel);
            app.ContrastUpperBoundLabel.HorizontalAlignment = 'center';
            app.ContrastUpperBoundLabel.WordWrap = 'on';
            app.ContrastUpperBoundLabel.Position = [26 64 86 29];
            app.ContrastUpperBoundLabel.Text = 'Contrast   Upper Bound';

            % Create ContrastUpperBoundSlider
            app.ContrastUpperBoundSlider = uislider(app.BrightnessandContrastPanel);
            app.ContrastUpperBoundSlider.Limits = [0 1];
            app.ContrastUpperBoundSlider.ValueChangingFcn = createCallbackFcn(app, @ContrastUpperBoundSliderValueChanging, true);
            app.ContrastUpperBoundSlider.Position = [130 83 150 3];
            app.ContrastUpperBoundSlider.Value = 1;

            % Create MaximizeContrastButton
            app.MaximizeContrastButton = uibutton(app.BrightnessandContrastPanel, 'push');
            app.MaximizeContrastButton.ButtonPushedFcn = createCallbackFcn(app, @MaximizeContrastButtonPushed, true);
            app.MaximizeContrastButton.Position = [139 21 115 22];
            app.MaximizeContrastButton.Text = 'Maximize Contrast';

            % Create ContrastLamp
            app.ContrastLamp = uilamp(app.BrightnessandContrastPanel);
            app.ContrastLamp.Enable = 'off';
            app.ContrastLamp.Position = [259 24 16 16];

            % Show the figure after all components are created
            app.AdjustImageUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = AdjustImageWindow_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.AdjustImageUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.AdjustImageUIFigure)
        end
    end
end