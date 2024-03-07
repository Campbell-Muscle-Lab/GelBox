classdef AdjustImageWindow_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        AdjustImageUIFigure             matlab.ui.Figure
        BrightnessandContrastPanel      matlab.ui.container.Panel
        MaximizeContrastCheckBox        matlab.ui.control.CheckBox
        MaximumInPixelIntensitySlider   matlab.ui.control.Slider
        MaximumInPixelIntensitySliderLabel  matlab.ui.control.Label
        BrightnessSlider                matlab.ui.control.Slider
        BrightnessSliderLabel           matlab.ui.control.Label
        MinimumInPixelIntensitySlider   matlab.ui.control.Slider
        MinimumInPixelIntensitySliderLabel  matlab.ui.control.Label
        adjusted_image_hist_working     matlab.ui.control.UIAxes
        in_out_axis                     matlab.ui.control.UIAxes
        AdjustedImagePanel              matlab.ui.container.Panel
        RevertChangesButton             matlab.ui.control.Button
        AcceptChangesButton             matlab.ui.control.Button
        adjusted_image_axis             matlab.ui.control.UIAxes
        OriginalImagePanel              matlab.ui.container.Panel
        RotateImageDegreesSpinner       matlab.ui.control.Spinner
        RotateImageDegreesSpinnerLabel  matlab.ui.control.Label
        InvertImageCheckBox             matlab.ui.control.CheckBox
        CropImageButton                 matlab.ui.control.Button
        original_image_axis             matlab.ui.control.UIAxes
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
        in_out % Description
        inverted_image % Description
        col = {'r','b'} % Description
        overlay_original_image % Description
        im_rotation_angle = [] % Description
    end

    methods (Access = public)

        function UpdateAdjustedImageDisplay(app,image)

            cla(app.adjusted_image_axis);
            cla(app.adjusted_image_hist_working);

            sat_mask = [];
            if app.inverted_image
                sat_mask = image == 0;
            else
                sat_mask = image == 1;
            end

            im_b = imoverlay(image,sat_mask,'red');
            center_image_with_preserved_aspect_ratio( ...
                im_b, ...
                app.adjusted_image_axis,[0 1]);
            colorbar(app.adjusted_image_axis)

            [counts, edges] = histcounts(image,linspace(0,1,101));
            yyaxis(app.adjusted_image_hist_working,'right')
            app.adjusted_image_hist_working.YAxis(2).Color = app.col{2};
            app.adjusted_image_hist_working.YAxis(1).Color = app.col{1};
            bar(app.adjusted_image_hist_working,...
                edges(1:(end-1)), counts,'EdgeColor',app.col{2},'FaceAlpha',0);
            xlim(app.adjusted_image_hist_working,[-0.01 1.01])
            OverlayOriginalHistogram(app)
        end

        function LoadImageAdjustments(app)

            app.original_image = app.GelBoxApp.gel_data.image.original_image;
            bit_d = class(app.original_image);

            app.bit_d_val = double(intmax(bit_d))+1;
            app.original_image = double(app.original_image)./(app.bit_d_val);

            app.overlay_original_image = app.original_image;
            center_image_with_preserved_aspect_ratio( ...
                app.original_image, ...
                app.original_image_axis,[]);

            col_bar = colorbar(app.original_image_axis);
            ax_lim(1) = col_bar.Limits(1) - 0.0001;
            ax_lim(2) = col_bar.Limits(2) + 0.0001;

            app.in_out = linspace(0,1);
            app.in_out = imadjust(app.in_out,[app.GelBoxApp.gel_data.settings.image_adjustments.contrast_lower...
                app.GelBoxApp.gel_data.settings.image_adjustments.contrast_upper],[]);
            plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
            OverlayOriginalHistogram(app)
            if app.GelBoxApp.gel_data.settings.image_adjustments.max_contrast
                app.MaximizeContrastCheckBox.Value = 1;
            end
            app.adjusted_image = app.GelBoxApp.gel_data.image.adjusted_image;
            bit_d = class(app.adjusted_image);

            app.bit_d_val = double(intmax(bit_d))+1;
            app.adjusted_image = double(app.adjusted_image)./(app.bit_d_val);


            app.BrightnessSlider.Value = app.GelBoxApp.gel_data.settings.image_adjustments.brightness;
            app.MinimumInPixelIntensitySlider.Value = app.GelBoxApp.gel_data.settings.image_adjustments.contrast_lower;
            app.MaximumInPixelIntensitySlider.Value = app.GelBoxApp.gel_data.settings.image_adjustments.contrast_upper;
            
            app.cropped_image = [];
            if ~isempty(app.GelBoxApp.gel_data.settings.image_adjustments.crop_pos)
                app.im_crop_box = images.roi.Rectangle(app.original_image_axis, ...
                    'Position',app.GelBoxApp.gel_data.settings.image_adjustments.crop_pos);
                app.im_crop_box.FaceAlpha = 0;
                app.im_crop_box.Color = [0 1 0];
                app.crop_pos = app.im_crop_box.Position;
                if isfield(app.GelBoxApp.gel_data.settings,'invert_status')
                    if app.GelBoxApp.gel_data.settings.invert_status
                        app.original_image = imcomplement(app.original_image);
                    end
                end
                app.cropped_image = imcrop(app.original_image,app.crop_pos);
                addlistener(app.im_crop_box,"MovingROI", ...
                    @(src,evt) crop_box_position_2(evt));
            end

            app.RevertChangesButton.Enable = 'off';

            %new status setting
            if ~isfield(app.GelBoxApp.gel_data.settings,'invert_status')
                app.InvertImageCheckBox.Value = 1;
                app.GelBoxApp.gel_data.settings.invert_status = 1;
                app.inverted_image = 1;
            elseif app.GelBoxApp.gel_data.settings.invert_status
                app.InvertImageCheckBox.Value = 1;
                app.inverted_image = 1;
                if isempty(app.cropped_image)
                    InvertImageCheckBoxValueChanged(app)
                end
            end
            
            try
                app.RotateImageDegreesSpinner.Value = ...
                    app.GelBoxApp.gel_data.settings.image_adjustments.im_rotation;
            catch
                app.RotateImageDegreesSpinner.Value = 0;
            end
           
            
            UpdateAdjustedImageDisplay(app,app.adjusted_image)
            
%             if app.GelBoxApp.loaded_analysis
%                 comp_names = {'CropImageButton','AcceptChangesButton',...
%                     'BrightnessSlider',...
%                     'MaximumInPixelIntensitySlider',...
%                     'MinimumInPixelIntensitySlider',...
%                     'MaximizeContrastButton',...
%                     'InvertImageCheckBox'};
%                 
%                 for i = 1 : numel(comp_names)
%                     app.(comp_names{i}).Enable = 'off';
%                 end
%                 app.im_crop_box.InteractionsAllowed = 'none';
%             end

            function crop_box_position_2(evt);
                app.cropped_image = [];
                app.crop_pos = app.im_crop_box.Position;
                app.cropped_image = imcrop(app.original_image,app.crop_pos);
%                 if app.inverted_image
%                     app.cropped_image = imcomplement(app.cropped_image);
%                 end
                app.BrightnessSlider.Value = 0;
                app.MinimumInPixelIntensitySlider.Value = 0;
                app.MaximumInPixelIntensitySlider.Value = 1;
                app.in_out = linspace(0,1);
                app.adjusted_image = [];
                app.adjusted_image_2 = [];
                plot(app.in_out_axis,app.in_out,app.in_out,'LineWidth',2,'Color','k')
                UpdateAdjustedImageDisplay(app,app.cropped_image);

            end

        end

        function OverlayOriginalHistogram(app)
            or_im = app.overlay_original_image;
            [counts, edges] = histcounts(or_im,linspace(0,1,101));
            yyaxis(app.adjusted_image_hist_working,'left')
            app.adjusted_image_hist_working.YAxis(1).Color = app.col{1};
            app.adjusted_image_hist_working.YAxis(2).Color = app.col{2};
            bar(app.adjusted_image_hist_working,...
                edges(1:(end-1)), counts,'EdgeColor',app.col{1},'FaceAlpha',0);
            xlim(app.adjusted_image_hist_working,[-0.01 1.01])
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, caller)
            movegui(app.AdjustImageUIFigure,'center')
            colormap(app.AdjustImageUIFigure, 'gray');
            xlim(app.adjusted_image_hist_working,[0 1]);
            xlim(app.adjusted_image_hist_working,[0 1]);

            app.GelBoxApp = caller;
            app.GelBoxApp.gel_data.invert_status = app.inverted_image;
            yyaxis(app.adjusted_image_hist_working,'left')
            ylabel(app.adjusted_image_hist_working,...
                {'Number of Pixels','(Original Image)'})
            yyaxis(app.adjusted_image_hist_working,'right')
            ylabel(app.adjusted_image_hist_working,...
                {'Number of Pixels','(Adjusted Image)'})

            app.adjusted_image_hist_working.YAxis(1).Exponent = 0;
            app.adjusted_image_hist_working.YAxis(2).Exponent = 0;
            %             app.hist_in.YAxis.Exponent = 0;
            %             app.hist_out.XAxis.Exponent = 0;

            if isfield(app.GelBoxApp.gel_data,'settings')
                LoadImageAdjustments(app)
            else
                app.original_image = caller.gel_data.image.im_data;

                bit_d = class(app.original_image);

                app.bit_d_val = double(intmax(bit_d))+1;
                app.original_image = double(app.original_image)./(app.bit_d_val);

                app.overlay_original_image = app.original_image;
                center_image_with_preserved_aspect_ratio( ...
                    app.original_image, ...
                    app.original_image_axis,[]);

                col_bar = colorbar(app.original_image_axis);
                ax_lim(1) = col_bar.Limits(1) - 0.0001;
                ax_lim(2) = col_bar.Limits(2) + 0.0001;
                app.in_out = linspace(0,1);
                plot(app.in_out_axis,app.in_out,app.in_out,'LineWidth',2,'Color','k')
                OverlayOriginalHistogram(app)
            end

        end

        % Button pushed function: CropImageButton
        function CropImageButtonPushed(app, event)
            app.im_crop_box = drawrectangle(app.original_image_axis);
            app.im_crop_box.FaceAlpha = 0;
            app.im_crop_box.Color = [0 1 0];
            app.crop_pos = app.im_crop_box.Position;
            if isempty(app.cropped_image) && ~isempty(app.adjusted_image)
                app.adjusted_image = [];
            end
            app.cropped_image = imcrop(app.original_image,app.crop_pos);
            addlistener(app.im_crop_box,"MovingROI", ...
                @(src,evt) crop_box_position(evt));
            

            UpdateAdjustedImageDisplay(app,app.cropped_image)

            function crop_box_position(evt);
                app.cropped_image = [];
                app.crop_pos = app.im_crop_box.Position;
                app.cropped_image = imcrop(app.overlay_original_image,app.crop_pos);
                if app.inverted_image
                    app.cropped_image = imcomplement(app.cropped_image);
                end
                app.BrightnessSlider.Value = 0;
                app.MinimumInPixelIntensitySlider.Value = 0;
                app.MaximumInPixelIntensitySlider.Value = 1;
                app.in_out = linspace(0,1);
                app.adjusted_image = [];
                app.adjusted_image_2 = [];
                plot(app.in_out_axis,app.in_out,app.in_out,'LineWidth',2,'Color','k')
                UpdateAdjustedImageDisplay(app,app.cropped_image);
            end

        end

        % Value changing function: BrightnessSlider
        function BrightnessSliderValueChanging(app, event)
            changingValue = event.Value;
            if app.MinimumInPixelIntensitySlider.Value == 0 && app.MaximumInPixelIntensitySlider.Value == 1
                if changingValue < 0
                    if isempty(app.cropped_image)
                        app.cropped_image = app.original_image;
                    end
                    app.adjusted_image = imadjust(app.cropped_image,[0 1], [0 1+changingValue]);
                    app.in_out = linspace(0,1);
                    app.in_out = imadjust(app.in_out,[0 1], [0 1+changingValue]);
                    cla(app.in_out_axis)
                    plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
                    xlim(app.in_out_axis,[0 1])
                    ylim(app.in_out_axis,[0 1])
                    app.MinimumInPixelIntensitySlider.Value = 0;
                    app.MaximumInPixelIntensitySlider.Value = 1+changingValue;
                    UpdateAdjustedImageDisplay(app,app.adjusted_image);
                else
                    if isempty(app.cropped_image)
                        app.cropped_image = app.original_image;
                    end
                    app.adjusted_image = imadjust(app.cropped_image,[0 1], [changingValue 1]);
                    app.in_out = linspace(0,1);
                    app.in_out = imadjust(app.in_out,[0 1], [changingValue 1]);
                    cla(app.in_out_axis)
                    plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
                    xlim(app.in_out_axis,[0 1])
                    ylim(app.in_out_axis,[0 1])
                    app.MinimumInPixelIntensitySlider.Value = changingValue;
                    app.MaximumInPixelIntensitySlider.Value = 1;
                    UpdateAdjustedImageDisplay(app,app.adjusted_image);
                end
            else
                if changingValue < 0
                    app.adjusted_image_2 = imadjust(app.adjusted_image,[0 1], [0 1+changingValue]);
                    app.in_out = linspace(0,1);
                    app.in_out = imadjust(app.in_out,[0 1], [0 1+changingValue]);
                    cla(app.in_out_axis)
                    plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
                    xlim(app.in_out_axis,[0 1])
                    ylim(app.in_out_axis,[0 1])
                    app.MinimumInPixelIntensitySlider.Value = 0;
                    app.MaximumInPixelIntensitySlider.Value = 1+changingValue;
                    UpdateAdjustedImageDisplay(app,app.adjusted_image_2);
                else
                    app.adjusted_image_2 = imadjust(app.adjusted_image,[0 1], [changingValue 1]);
                    app.in_out = linspace(0,1);
                    app.in_out = imadjust(app.in_out,[0 1], [changingValue 1]);
                    cla(app.in_out_axis)
                    plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
                    xlim(app.in_out_axis,[0 1])
                    ylim(app.in_out_axis,[0 1])
                    app.MinimumInPixelIntensitySlider.Value = changingValue;
                    app.MaximumInPixelIntensitySlider.Value = 1;
                    UpdateAdjustedImageDisplay(app,app.adjusted_image_2);
                end
            end
        end

        % Value changing function: MinimumInPixelIntensitySlider
        function MinimumInPixelIntensitySliderValueChanging(app, event)
            changingValue = event.Value;
            if app.BrightnessSlider.Value == 0 && isempty(app.adjusted_image)
                if isempty(app.cropped_image)
                    im = app.original_image;
                else
                    im = app.cropped_image;
                end
                app.adjusted_image = imadjust(im,[changingValue app.MaximumInPixelIntensitySlider.Value],[]);
                app.in_out = linspace(0,1);
                app.in_out = imadjust(app.in_out,[changingValue app.MaximumInPixelIntensitySlider.Value],[]);
                cla(app.in_out_axis)
                plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
                xlim(app.in_out_axis,[0 1])
                ylim(app.in_out_axis,[0 1])
                UpdateAdjustedImageDisplay(app,app.adjusted_image);
            else
                app.adjusted_image_2 = imadjust(app.adjusted_image,[changingValue app.MaximumInPixelIntensitySlider.Value],[]);
                app.in_out = linspace(0,1);
                app.in_out = imadjust(app.in_out,[changingValue app.MaximumInPixelIntensitySlider.Value],[]);
                cla(app.in_out_axis)
                plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
                xlim(app.in_out_axis,[0 1])
                ylim(app.in_out_axis,[0 1])
                UpdateAdjustedImageDisplay(app,app.adjusted_image_2);
            end
        end

        % Value changing function: MaximumInPixelIntensitySlider
        function MaximumInPixelIntensitySliderValueChanging(app, event)
            changingValue = event.Value;
            if app.BrightnessSlider.Value == 0 && isempty(app.adjusted_image)
                if isempty(app.cropped_image)
                    im = app.original_image;
                else
                    im = app.cropped_image;
                end
                app.adjusted_image = imadjust(im,[app.MinimumInPixelIntensitySlider.Value changingValue],[]);
                app.in_out = linspace(0,1);
                app.in_out = imadjust(app.in_out,[app.MinimumInPixelIntensitySlider.Value changingValue],[]);
                cla(app.in_out_axis)
                plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
                xlim(app.in_out_axis,[0 1])
                ylim(app.in_out_axis,[0 1])
                UpdateAdjustedImageDisplay(app,app.adjusted_image);
            else
                app.adjusted_image_2 = imadjust(app.adjusted_image,[app.MinimumInPixelIntensitySlider.Value changingValue],[]);
                app.in_out = linspace(0,1);
                app.in_out = imadjust(app.in_out,[app.MinimumInPixelIntensitySlider.Value changingValue],[]);
                cla(app.in_out_axis)
                plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
                xlim(app.in_out_axis,[0 1])
                ylim(app.in_out_axis,[0 1])
                UpdateAdjustedImageDisplay(app,app.adjusted_image_2);
            end
        end

        % Button pushed function: RevertChangesButton
        function RevertChangesButtonPushed(app, event)
            app.im_rotation_angle = [];
            app.cropped_image = [];
            app.adjusted_image = [];
            app.adjusted_image_2 = [];
            app.RotateImageDegreesSpinner.Value = 0;
            app.MinimumInPixelIntensitySlider.Value = 0;
            app.MaximumInPixelIntensitySlider.Value = 1;
            app.BrightnessSlider.Value = 0;
            app.in_out = linspace(0,1);
            app.MaximizeContrastCheckBox.Value = 0;
            app.max_contrast = false;
            plot(app.in_out_axis,app.in_out,app.in_out,'LineWidth',2,'Color','k')
            if app.inverted_image
                app.original_image = app.overlay_original_image;
                app.inverted_image = 0;
                app.InvertImageCheckBox.Value = 0;
            end
            delete(app.im_crop_box)
            UpdateAdjustedImageDisplay(app,app.original_image)

        end

        % Button pushed function: AcceptChangesButton
        function AcceptChangesButtonPushed(app, event)

            app.GelBoxApp.gel_data.settings.image_adjustments.crop_pos = app.crop_pos;
            app.GelBoxApp.gel_data.settings.image_adjustments.brightness = app.BrightnessSlider.Value;
            app.GelBoxApp.gel_data.settings.image_adjustments.contrast_lower = app.MinimumInPixelIntensitySlider.Value;
            app.GelBoxApp.gel_data.settings.image_adjustments.contrast_upper = app.MaximumInPixelIntensitySlider.Value;
            app.GelBoxApp.gel_data.settings.image_adjustments.max_contrast = app.max_contrast;
            app.GelBoxApp.gel_data.settings.image_adjustments.im_rotation = app.RotateImageDegreesSpinner.Value;
            app.GelBoxApp.gel_data.settings.invert_status = app.inverted_image;
            if isempty(app.cropped_image) && isempty(app.im_rotation_angle)
                if app.inverted_image
                    app.GelBoxApp.gel_data.image.adjusted_image = imcomplement(app.GelBoxApp.gel_data.image.original_image);
                else
                    app.GelBoxApp.gel_data.image.adjusted_image = [];
                end
            elseif isempty(app.adjusted_image)
                app.GelBoxApp.gel_data.image.adjusted_image = app.cropped_image.*app.bit_d_val;
            elseif isempty(app.adjusted_image_2)
                app.GelBoxApp.gel_data.image.adjusted_image = app.adjusted_image.*app.bit_d_val;
            else
                app.GelBoxApp.gel_data.image.adjusted_image = app.adjusted_image_2.*app.bit_d_val;
            end
            
            app.GelBoxApp.gel_data.image.adjusted_image = cast(app.GelBoxApp.gel_data.image.adjusted_image,...
                class(app.GelBoxApp.gel_data.image.im_data));

            ImageAdjusted(app.GelBoxApp)
            delete(app)
        end

        % Value changed function: InvertImageCheckBox
        function InvertImageCheckBoxValueChanged(app, event)

            value = app.InvertImageCheckBox.Value;
            if isempty(app.cropped_image) && isempty(app.im_rotation_angle)
                app.original_image = imcomplement(app.original_image);
                im = app.original_image;
            elseif isempty(app.adjusted_image)
                app.cropped_image = imcomplement(app.cropped_image);
                im = app.cropped_image;
            elseif isempty(app.adjusted_image_2)
                app.adjusted_image = imcomplement(app.adjusted_image);
                im = app.adjusted_image;
            else
                app.adjusted_image_2 = imcomplement(app.adjusted_image_2);
                im = app.adjusted_image_2;
            end
            app.inverted_image = value;
            UpdateAdjustedImageDisplay(app,im)

        end

        % Value changed function: MaximizeContrastCheckBox
        function MaximizeContrastCheckBoxValueChanged(app, event)
            value = app.MaximizeContrastCheckBox.Value;
            if value
                if isempty(app.adjusted_image)
                    if isempty(app.cropped_image)
                        app.cropped_image = app.original_image;
                    end
                    app.adjusted_image = imadjust(app.cropped_image,stretchlim(app.cropped_image,[0 1]),[0 1]);
                    app.in_out = imadjust(app.in_out,stretchlim(app.in_out,[0 1]),[0 1]);
                    cla(app.in_out_axis)
                    plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
                    xlim(app.in_out_axis,[0 1])
                    ylim(app.in_out_axis,[0 1])
                    UpdateAdjustedImageDisplay(app,app.adjusted_image);
                else
                    app.adjusted_image_2 = imadjust(app.adjusted_image,stretchlim(app.adjusted_image,[0 1]),[0 1]);
                    app.in_out = imadjust(app.in_out,stretchlim(app.in_out,[0 1]),[0 1]);
                    cla(app.in_out_axis)
                    plot(app.in_out_axis,linspace(0,1),app.in_out,'LineWidth',2,'Color','k')
                    xlim(app.in_out_axis,[0 1])
                    ylim(app.in_out_axis,[0 1])
                    UpdateAdjustedImageDisplay(app,app.adjusted_image_2);
                end
            else
                if ~isempty(app.adjusted_image_2)
                    app.adjusted_image_2 = [];
                    UpdateAdjustedImageDisplay(app,app.adjusted_image);
                elseif ~isempty(app.adjusted_image)
                    app.adjusted_image = [];
                    if isempty(app.cropped_image)
                        app.cropped_image = app.original_image;
                    end
                    UpdateAdjustedImageDisplay(app,app.cropped_image);
                end
            end

            app.max_contrast = true;
        end

        % Value changing function: RotateImageDegreesSpinner
        function RotateImageDegreesSpinnerValueChanging(app, event)
            changingValue = event.Value;
            if app.im_rotation_angle
                rot = changingValue - app.im_rotation_angle;
            else
                rot = changingValue;
            end
            if app.inverted_image
                pad_val = 1;
            else
                pad_val = 0;
            end
            rot_mask = [];

            if isempty(app.adjusted_image)
                if isempty(app.cropped_image)
                    app.adjusted_image = imrotate(app.original_image,rot,'crop');
                    rot_mask = ~imrotate(app.original_image,rot,'crop');
                else
                    app.adjusted_image = imrotate(app.cropped_image,rot,'crop');
                    rot_mask = ~imrotate(app.cropped_image,rot,'crop');
                end

                app.adjusted_image(rot_mask) = pad_val;
                im = app.adjusted_image;

                app.im_rotation_angle = changingValue;

            elseif isempty(app.adjusted_image_2)
                app.adjusted_image_2 = imrotate(app.adjusted_image,rot,'crop');
                rot_mask = ~imrotate(app.adjusted_image,rot,'crop');
                app.adjusted_image_2(rot_mask) = pad_val;
                im = app.adjusted_image_2;

                app.im_rotation_angle = changingValue;
            else
                app.adjusted_image_2 = imrotate(app.adjusted_image_2,rot,'crop');
                rot_mask = ~imrotate(app.adjusted_image_2,rot,'crop');
                app.adjusted_image_2(rot_mask) = pad_val;
                im = app.adjusted_image_2;


                app.im_rotation_angle = changingValue;
            end
            UpdateAdjustedImageDisplay(app,im);


        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create AdjustImageUIFigure and hide until all components are created
            app.AdjustImageUIFigure = uifigure('Visible', 'off');
            app.AdjustImageUIFigure.Position = [100 100 1690 582];
            app.AdjustImageUIFigure.Name = 'Adjust Image';

            % Create OriginalImagePanel
            app.OriginalImagePanel = uipanel(app.AdjustImageUIFigure);
            app.OriginalImagePanel.Title = 'Original Image';
            app.OriginalImagePanel.Position = [8 18 558 553];

            % Create original_image_axis
            app.original_image_axis = uiaxes(app.OriginalImagePanel);
            app.original_image_axis.XTick = [];
            app.original_image_axis.YTick = [];
            app.original_image_axis.Box = 'on';
            app.original_image_axis.Position = [8 0 540 484];

            % Create CropImageButton
            app.CropImageButton = uibutton(app.OriginalImagePanel, 'push');
            app.CropImageButton.ButtonPushedFcn = createCallbackFcn(app, @CropImageButtonPushed, true);
            app.CropImageButton.Position = [12 495 100 22];
            app.CropImageButton.Text = 'Crop Image';

            % Create InvertImageCheckBox
            app.InvertImageCheckBox = uicheckbox(app.OriginalImagePanel);
            app.InvertImageCheckBox.ValueChangedFcn = createCallbackFcn(app, @InvertImageCheckBoxValueChanged, true);
            app.InvertImageCheckBox.Text = 'Invert Image';
            app.InvertImageCheckBox.Position = [126 495 89 22];

            % Create RotateImageDegreesSpinnerLabel
            app.RotateImageDegreesSpinnerLabel = uilabel(app.OriginalImagePanel);
            app.RotateImageDegreesSpinnerLabel.HorizontalAlignment = 'right';
            app.RotateImageDegreesSpinnerLabel.Position = [222 495 134 22];
            app.RotateImageDegreesSpinnerLabel.Text = 'Rotate Image (Degrees)';

            % Create RotateImageDegreesSpinner
            app.RotateImageDegreesSpinner = uispinner(app.OriginalImagePanel);
            app.RotateImageDegreesSpinner.ValueChangingFcn = createCallbackFcn(app, @RotateImageDegreesSpinnerValueChanging, true);
            app.RotateImageDegreesSpinner.Position = [363 495 50 22];

            % Create AdjustedImagePanel
            app.AdjustedImagePanel = uipanel(app.AdjustImageUIFigure);
            app.AdjustedImagePanel.Title = 'Adjusted Image';
            app.AdjustedImagePanel.Position = [1105 18 558 553];

            % Create adjusted_image_axis
            app.adjusted_image_axis = uiaxes(app.AdjustedImagePanel);
            app.adjusted_image_axis.XTick = [];
            app.adjusted_image_axis.YTick = [];
            app.adjusted_image_axis.Box = 'on';
            app.adjusted_image_axis.Position = [10 7 540 478];

            % Create AcceptChangesButton
            app.AcceptChangesButton = uibutton(app.AdjustedImagePanel, 'push');
            app.AcceptChangesButton.ButtonPushedFcn = createCallbackFcn(app, @AcceptChangesButtonPushed, true);
            app.AcceptChangesButton.Position = [10 495 104 22];
            app.AcceptChangesButton.Text = 'Accept Changes';

            % Create RevertChangesButton
            app.RevertChangesButton = uibutton(app.AdjustedImagePanel, 'push');
            app.RevertChangesButton.ButtonPushedFcn = createCallbackFcn(app, @RevertChangesButtonPushed, true);
            app.RevertChangesButton.Position = [130 495 102 22];
            app.RevertChangesButton.Text = 'Revert Changes';

            % Create BrightnessandContrastPanel
            app.BrightnessandContrastPanel = uipanel(app.AdjustImageUIFigure);
            app.BrightnessandContrastPanel.Title = 'Brightness and Contrast';
            app.BrightnessandContrastPanel.Position = [577 18 519 552];

            % Create in_out_axis
            app.in_out_axis = uiaxes(app.BrightnessandContrastPanel);
            xlabel(app.in_out_axis, 'In (Pixel Intensity)')
            ylabel(app.in_out_axis, 'Out (Pixel Intensity)')
            zlabel(app.in_out_axis, 'Z')
            app.in_out_axis.Box = 'on';
            app.in_out_axis.Position = [222 27 235 235];

            % Create adjusted_image_hist_working
            app.adjusted_image_hist_working = uiaxes(app.BrightnessandContrastPanel);
            title(app.adjusted_image_hist_working, 'Intensity Histogram')
            xlabel(app.adjusted_image_hist_working, 'Pixel Intensity')
            app.adjusted_image_hist_working.Box = 'on';
            app.adjusted_image_hist_working.Position = [32 269 458 256];

            % Create MinimumInPixelIntensitySliderLabel
            app.MinimumInPixelIntensitySliderLabel = uilabel(app.BrightnessandContrastPanel);
            app.MinimumInPixelIntensitySliderLabel.HorizontalAlignment = 'center';
            app.MinimumInPixelIntensitySliderLabel.WordWrap = 'on';
            app.MinimumInPixelIntensitySliderLabel.Position = [30 163 156 29];
            app.MinimumInPixelIntensitySliderLabel.Text = 'Minimum In (Pixel Intensity)';

            % Create MinimumInPixelIntensitySlider
            app.MinimumInPixelIntensitySlider = uislider(app.BrightnessandContrastPanel);
            app.MinimumInPixelIntensitySlider.Limits = [0 1];
            app.MinimumInPixelIntensitySlider.ValueChangingFcn = createCallbackFcn(app, @MinimumInPixelIntensitySliderValueChanging, true);
            app.MinimumInPixelIntensitySlider.Position = [32 151 150 3];

            % Create BrightnessSliderLabel
            app.BrightnessSliderLabel = uilabel(app.BrightnessandContrastPanel);
            app.BrightnessSliderLabel.HorizontalAlignment = 'center';
            app.BrightnessSliderLabel.WordWrap = 'on';
            app.BrightnessSliderLabel.Position = [77 235 57 30];
            app.BrightnessSliderLabel.Text = 'Brightness';

            % Create BrightnessSlider
            app.BrightnessSlider = uislider(app.BrightnessandContrastPanel);
            app.BrightnessSlider.Limits = [-1 1];
            app.BrightnessSlider.MajorTicks = [-1 -0.6 -0.2 0 0.2 0.6 1];
            app.BrightnessSlider.MajorTickLabels = {'Darken', '', '', '', '', '', 'Brighten'};
            app.BrightnessSlider.ValueChangingFcn = createCallbackFcn(app, @BrightnessSliderValueChanging, true);
            app.BrightnessSlider.Position = [29 226 153 3];

            % Create MaximumInPixelIntensitySliderLabel
            app.MaximumInPixelIntensitySliderLabel = uilabel(app.BrightnessandContrastPanel);
            app.MaximumInPixelIntensitySliderLabel.HorizontalAlignment = 'center';
            app.MaximumInPixelIntensitySliderLabel.WordWrap = 'on';
            app.MaximumInPixelIntensitySliderLabel.Position = [29 80 161 29];
            app.MaximumInPixelIntensitySliderLabel.Text = 'Maximum In (Pixel Intensity)';

            % Create MaximumInPixelIntensitySlider
            app.MaximumInPixelIntensitySlider = uislider(app.BrightnessandContrastPanel);
            app.MaximumInPixelIntensitySlider.Limits = [0 1];
            app.MaximumInPixelIntensitySlider.ValueChangingFcn = createCallbackFcn(app, @MaximumInPixelIntensitySliderValueChanging, true);
            app.MaximumInPixelIntensitySlider.Position = [28 73 156 3];
            app.MaximumInPixelIntensitySlider.Value = 1;

            % Create MaximizeContrastCheckBox
            app.MaximizeContrastCheckBox = uicheckbox(app.BrightnessandContrastPanel);
            app.MaximizeContrastCheckBox.ValueChangedFcn = createCallbackFcn(app, @MaximizeContrastCheckBoxValueChanged, true);
            app.MaximizeContrastCheckBox.Text = 'Maximize Contrast';
            app.MaximizeContrastCheckBox.Position = [51 7 121 22];

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