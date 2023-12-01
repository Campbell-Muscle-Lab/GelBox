classdef SummaryPlotWindow_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        SummaryPlotUIFigure           matlab.ui.Figure
        ControlLabelsCommaSeparatedEditField  matlab.ui.control.EditField
        ControlLabelsCommaSeparatedEditFieldLabel  matlab.ui.control.Label
        BandLabelsPanel               matlab.ui.container.Panel
        LaneNumbersCommaSeparatedEditField  matlab.ui.control.EditField
        LaneNumbersCommaSeparatedEditFieldLabel  matlab.ui.control.Label
        ControlLanesCheckBox          matlab.ui.control.CheckBox
        BandLabelsCommaSeparetedEditField  matlab.ui.control.EditField
        BandLabelsCommaSeparetedEditFieldLabel  matlab.ui.control.Label
        FileOptionsPanel              matlab.ui.container.Panel
        FigureFileNameField           matlab.ui.control.EditField
        FigureFileNameEditFieldLabel  matlab.ui.control.Label
        OutputPathField               matlab.ui.control.EditField
        SelectOutputFolderButton      matlab.ui.control.Button
        GenerateSummaryPlotsButton    matlab.ui.control.Button
    end


    properties (Access = private)
        GelBoxApp % Description
    end

    properties (Access = public)
        output_path % Description
        box_layout_fname
        f_title
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, caller)

            movegui(app.SummaryPlotUIFigure,'center')


            app.GelBoxApp = caller;


        end

        % Button pushed function: SelectOutputFolderButton
        function SelectOutputFolderButtonPushed(app, event)
            app.output_path = uigetdir2('','Select Output Folder');
            if (app.output_path~=0)
            app.OutputPathField.Value = app.output_path;
            end
        end

        % Value changed function: FigureFileNameField
        function FigureFileNameFieldValueChanged(app, event)
            fname = app.FigureFileNameField.Value;
            app.f_title = fname;
            app.box_layout_fname = sprintf('%s\\%s_box_layout', ...
                app.output_path,fname);
            app.BandLabelsCommaSeparetedEditField.Enable = 1;
            app.BandLabelsCommaSeparetedEditFieldLabel.Enable = 1;
            app.ControlLanesCheckBox.Enable = 1;

        end

        % Button pushed function: GenerateSummaryPlotsButton
        function GenerateSummaryPlotsButtonPushed(app, event)
            
            user_label_entry = app.BandLabelsCommaSeparetedEditField.Value;
            
            if app.ControlLanesCheckBox.Value
                
                control_lane_entry = app.LaneNumbersCommaSeparatedEditField.Value;
                
                control_lanes = split(control_lane_entry,",");
                
                control_lanes = str2double(control_lanes(:,1));
                
                control_label_entry = app.ControlLabelsCommaSeparatedEditField.Value;
                
                control_labels = strsplit(control_label_entry);
                
                for i = 1 : numel(control_labels)
                control_labels{i} = regexprep(control_labels{i},',','');
                control_labels{i} = regexprep(control_labels{i},' ','');
                end
            else
                control_lanes = [];
            end
                            
            sample_labels = split(user_label_entry,",");
            
            for i = 1 : numel(sample_labels)
                sample_labels{i} = regexprep(sample_labels{i},',','');
                sample_labels{i} = regexprep(sample_labels{i},' ','');
            end
            
            num_of_bands = numel(sample_labels);
            

            max_x = [];
            max_x_fit = [];
            scaling = {'same_scale','auto_scale'};
            summary = app.GelBoxApp.gel_data.summary;
            [m,n] = size(summary);
            for i = 1 : n
                max_x(i) = max(summary(i).x);
                max_x_fit(i) = max(summary(i).x_fit);
                min_x_fit(i) = min(summary(i).x_fit);
            end

            max_x = ceil(max(max_x)/100)*100;
            max_x_fit = ceil(max(max_x_fit)/100)*100;
            min_x_fit = ceil(min(min_x_fit)/-100)*-100;


            box_positions = [];

            for i = 1 : length(summary)
                box_positions(i,:) = summary(i).box_position;
            end



            im_a = app.GelBoxApp.gel_data.image.im_data;
            im_or = app.GelBoxApp.gel_data.image.original_image;


            if n == 1
                no_of_panels_wide = 3;
            else
                no_of_panels_wide = 6;
            end

            total_figures = n*3;
            fig_handles = 1:total_figures;
            no_of_panels_high = 2 + (total_figures + ...
                mod(total_figures,no_of_panels_wide))/no_of_panels_wide;

            left_pads = 0.05 * ones(1,no_of_panels_wide);
            left_pads([1]) = 0.4;
            left_pads([2 5]) = 0.0005;
            left_pads([4]) = 0.7;
            right_pads = 0.25 * ones(1,no_of_panels_wide);
            right_pads(3) = 0.5; 
            if no_of_panels_wide == 3
                right_pads(no_of_panels_wide) = 0.95;
            else
                right_pads(no_of_panels_wide) = 0.95;
            end

            omit_panels = [2:no_of_panels_wide];
            omit_panels = [omit_panels, no_of_panels_wide+2:2*no_of_panels_wide];

            if mod(total_figures,no_of_panels_wide) ~= 0
                res_total_figures = no_of_panels_high * no_of_panels_wide;
                res_total_figures = 1 : res_total_figures;
                omit_panels = [omit_panels ...
                    res_total_figures(end-2) ...
                    res_total_figures(end-1) ...
                    res_total_figures(end)];
            else
                omit_panels = omit_panels;
            end

            im_handles = 2*no_of_panels_wide + 1:3:total_figures + 2*no_of_panels_wide;
            back_fit_handles = 2*no_of_panels_wide + 2:3:total_figures + 2*no_of_panels_wide;
            fit_handles = 2*no_of_panels_wide + 3:3:total_figures + 2*no_of_panels_wide;

            right_subplot_adjustments = zeros(1, no_of_panels_wide*no_of_panels_high);
            right_subplot_adjustments(1) = -no_of_panels_wide + 0.005;
            right_subplot_adjustments(1+no_of_panels_wide) = -no_of_panels_wide + 0.005;
            relative_row_heights = ones(1,no_of_panels_high);
            relative_row_heights(1) = 3;
            relative_row_heights(2) = 1;
            
            for fig_no = 1 : 2
                colormap('gray')
                sp = initialise_publication_quality_figure( ...
                    'no_of_panels_wide', no_of_panels_wide, ...
                    'no_of_panels_high', no_of_panels_high, ...
                    'top_margin', 0.5, ...
                    'bottom_margin', 0.15, ...
                    'right_margin', 0.15, ...
                    'individual_padding', 1, ...
                    'left_pads', repmat(left_pads,[1 no_of_panels_high]), ...
                    'right_pads', repmat(right_pads,[1 no_of_panels_high]), ...
                    'axes_padding_top', 0.1, ...
                    'axes_padding_bottom',0.35, ...
                    'panel_label_font_size', 0, ...
                    'omit_panels', omit_panels, ...
                    'right_subplot_adjustments', right_subplot_adjustments,...
                    'relative_row_heights', relative_row_heights,...
                    'figure_handle',1);
                                
                subplot(sp(1))
                title('Original Image','FontSize',11)
                hold on
                center_image_with_preserved_aspect_ratio( ...
                    im_or, ...
                    sp(1),[])
                [m,n] = size(im_or);
                rectangle('Position',[1 1 n-1 m-1 ],'EdgeColor','k','LineWidth',1)
                
                subplot(sp(7))
                title('Adjusted Image','FontSize',11)
                hold on
                center_image_with_preserved_aspect_ratio( ...
                    im_a, ...
                    sp(7),[])

                for i = 1 : length(summary)
                    rectangle('Position',box_positions(i,:),...
                        'EdgeColor','g',...
                        'LineWidth',2)
                    text(box_positions(i,1)+box_positions(i,3)+17,box_positions(i,2)-50, ...
                        sprintf('%.0f',i),'FontWeight',"bold","FontSize",8);
                end

                [m,n] = size(im_a);
                rectangle('Position',[1 1 n-1 m-1 ],'EdgeColor','k','LineWidth',1)


                for i = 1 : length (im_handles)
                    h = subplot(sp(im_handles(i)));
                    colormap(h,"gray");
                    t1 = sprintf('Box %i',i);
                    title(t1, 'FontSize',7);
                    center_image_with_preserved_aspect_ratio( ...
                        summary(i).summary_inset,sp(im_handles(i)),[]);
                    x = size(summary(i).summary_inset,2) - size(summary(i).inset,2);
                    x = 0.5*x + 1;
                    y = size(summary(i).summary_inset,1) - size(summary(i).inset,1);
                    y = 0.5*y + 1;
                    rectangle('Position',[11 11 size(summary(i).inset,2) size(summary(i).inset,1)],'EdgeColor','g')
                    subplot(sp(back_fit_handles(i)))
                    plot(summary(i).x,summary(i).y,'k', 'LineWidth',1);

                    hold on
                    subplot(sp(fit_handles(i)))
                    plot(summary(i).x-(summary(i).x_back), ...
                        summary(i).y,'-.k', 'LineWidth',1);
                    hold on
                    
                    cont = 1;
                    
                    for j = 1 : size(summary(i).band,2)
                        
                        if ~isempty(control_lanes) && control_lanes(cont) == i
                            colors = pink(numel(control_lanes));
                            labels = control_labels{cont};
                            labels = string(labels);
                            cont = cont + 1;
                        else
                            colors = parula(numel(summary(i).area));
                            labels = sample_labels;
                            labels = string(labels);
                        end
                        
                        if numel(summary(i).area) < numel(labels)
                            for c = 2:numel(labels)
                            labels{c-1} = sample_labels{c};
                            end
                            labels = string(labels);
                        end
                        
                        subplot(sp(fit_handles(i)))
                        f(j) = patch(summary(i).band(:,j), ...
                            summary(i).y,colors(j,:),'FaceAlpha',0.35, ...
                            'EdgeColor',colors(j,:),'EdgeAlpha',0.05);
                        
                        rel_area(j) = (summary(i).area(j))/sum(summary(i).area);

                        subplot(sp(back_fit_handles(i)))
                        hold on

                        str_labels{j} = sprintf('%s = %.2f', ...
                            labels(j),rel_area(j));
                    end   

                    leg_labels = str_labels;
                    figs = f;

                    f_leg = sprintf('R^2 = %.3f',summary(i).r_squared);
                    f_color = '#1aff00';

                    subplot(sp(back_fit_handles(i)))
                    hold on
                    plot(summary(i).x_back,summary(i).y, ...
                        '-.m','LineWidth',1);
                    plot(summary(i).x_back + summary(i).x_fit, ...
                        summary(i).y,':',"LineWidth",1,'Color',f_color)

                    subplot(sp(fit_handles(i)))
                    hold on
                    plot(zeros(1,numel(summary(i).y)),summary(i).y, ...
                        '-.m','LineWidth',1);
                    f = plot(summary(i).x_fit, ...
                        summary(i).y,':',"LineWidth",1,'Color',f_color);
                    leg_labels{end+1} = f_leg;
                    figs(end+1) = f;
                    
                    back_method = regexp(app.GelBoxApp.gel_data.settings.background.method{i},'[^()]*','match');
                    if strcmp(back_method{2},'CSS')
                        fraction = app.GelBoxApp.gel_data.settings.background.css_fraction(i);
                        fraction_ix = ceil(fraction*numel(summary(i).x));
                        ix_1 = 1:fraction_ix;
                        ix_2 = numel(summary(i).x):-1:(numel(summary(i).x) - fraction_ix + 1);
                        c_x = [];
                        c_y = [];
                        c_x = [ix_1 ix_2];
                        c_y = summary(i).x(c_x);
                        s = app.GelBoxApp.gel_data.settings.background.css_smoothing(i);
                        n=0;
                        while (floor(s*10^n)~=s*10^n)
                            n=n+1;
                        end
                        formatting = sprintf('%s%i%s','%s: Fr (%%) = %.0f S = %.',n,'f');
                        figs(end+1) = plot(sp(back_fit_handles(i)),c_y,c_x,'ro','LineWidth',1,'MarkerSize',1);
                        leg_labels{end+1} = sprintf(formatting,back_method{2},fraction*100,s);
                    end

                    legendflex(figs, leg_labels, ...
                        'xscale',0.15, ...
                        'anchor',{'ne','ne'}, ...
                        'buffer',[90 0], ...
                        'padding',[1 1 2], ...
                        'FontSize',4.5, ...
                        'text_y_padding', 0);
                    rel_area = [];
                    str_labels = {};
                end
%                 sgtitle(app.f_title,'FontSize',15);                
                for i = 1 : length(im_handles)
                    improve_axes('axis_handle',sp(im_handles(i)),...
                        'y_axis_label','Pixels', ...
                        'y_tick_decimal_places',0,...
                        'y_label_offset',-0.15,...
                        'x_tick_decimal_places',0, ...
                        'y_tick_decimal_places',0, ...
                        'x_axis_label','Pixels', ...
                        'x_label_offset',-0.35,...
                        'x_tick_labels',{'1',string(size(summary(i).inset,2))},...
                        'x_tick_label_positions',[11 10+size(summary(i).inset,2)],...
                        'y_axis_offset',-0.05,'y_ticks',[size(summary(i).summary_inset,1) 1],...
                        'y_tick_labels',{'1',string(size(summary(i).inset,1))},...
                        'y_tick_label_positions',[size(summary(i).inset,1)+10 11],...
                        'clip_x_axis',1,...
                        'clip_y_axis',1,...
                        'tick_font_size',8,...
                        'label_font_size',7,...
                        'y_label_rotation',90)
                end

                for i = 1:length(fit_handles)

                    if fig_no == 1
                        x_ticks = [min_x_fit max_x_fit];
                    else
                        x_pow = ceil(log10(max(summary(i).x-summary(i).x_back)));
                        x_tick_rounder = 10^(x_pow - 1);
                        x_t_end = ceil(max(summary(i).x-summary(i).x_back) ...
                            /x_tick_rounder)*x_tick_rounder;
                        if min(summary(i).x-summary(i).x_back) < 0
                            x_tick_rounder = -10^(x_pow - 2)*0.25;
                        else
                            x_tick_rounder = 10^(x_pow - 2)*0.25;
                        end
                        x_t_beginning = ceil(min(summary(i).x-summary(i).x_back)/x_tick_rounder)*x_tick_rounder;
                        x_ticks = [x_t_beginning x_t_end];
                    end

                    improve_axes('axis_handle',sp(fit_handles(i)),...
                        'x_tick_decimal_places',0, ...
                        'y_tick_decimal_places',0,...
                        'y_axis_label',{'Pixx_ticks = [min_x_fit max_x_fit];x_ticks = [min_x_fit max_x_fit];els'},'y_label_offset',-0.4,...
                        'x_label_offset',-0.30,...
                        'tick_font_size',8,...
                        'label_font_size',8,...
                        'x_axis_label',{'Opt.','Density (A.U.)'},...
                        'x_ticks', x_ticks,...
                        'y_ticks',[1 length(summary(i).x)],...
                        'title_font_size', 8, ...
                        'title_y_offset',1.2,...
                        'y_axis_off',1)
                end

                for i = 1:length(back_fit_handles)

                    if fig_no == 1
                        x_ticks = [0 max_x];
                    else
                        x_pow = ceil(log10(max(summary(i).x)));
                        x_tick_rounder = 10^(x_pow - 1);
                        x_t_end = ceil(max(summary(i).x)/x_tick_rounder)*x_tick_rounder;
                        x_ticks = [0 x_t_end];
                    end


                    improve_axes('axis_handle',sp(back_fit_handles(i)),...
                        'x_tick_decimal_places',0, ...
                        'y_tick_decimal_places',0,...
                        'y_axis_off',0,...
                        'x_label_offset',-0.30,...
                        'y_label_offset',-0.15,...
                        'tick_font_size',8,...
                        'label_font_size',8,...
                        'y_axis_label',{'Pixels'},...
                        'x_axis_label',{'Opt.','Density (A.U.)'},...
                        'x_ticks',x_ticks,...
                        'y_ticks',[1 length(summary(i).x)],...
                        'y_label_rotation',90)
                end
                fname = sprintf('%s_%s',app.box_layout_fname,scaling{fig_no});
                figure_export('output_file_string', fname, ...
                    'output_type', 'png');
            end
            delete(app)





        end

        % Value changed function: ControlLanesCheckBox
        function ControlLanesCheckBoxValueChanged(app, event)
            value = app.ControlLanesCheckBox.Value;
            if value == 1
                app.ControlLabelsCommaSeparatedEditField.Enable = 1;
                app.ControlLabelsCommaSeparatedEditFieldLabel.Enable = 1;
                app.LaneNumbersCommaSeparatedEditField.Enable = 1;
                app.LaneNumbersCommaSeparatedEditFieldLabel.Enable = 1;
            else
                app.ControlLabelsCommaSeparatedEditField.Enable = 0;
                app.ControlLabelsCommaSeparatedEditFieldLabel.Enable = 0;
                app.LaneNumbersCommaSeparatedEditField.Enable = 0;
                app.LaneNumbersCommaSeparatedEditFieldLabel.Enable = 0;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create SummaryPlotUIFigure and hide until all components are created
            app.SummaryPlotUIFigure = uifigure('Visible', 'off');
            app.SummaryPlotUIFigure.Position = [100 100 345 416];
            app.SummaryPlotUIFigure.Name = 'Summary Plot';

            % Create GenerateSummaryPlotsButton
            app.GenerateSummaryPlotsButton = uibutton(app.SummaryPlotUIFigure, 'push');
            app.GenerateSummaryPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateSummaryPlotsButtonPushed, true);
            app.GenerateSummaryPlotsButton.Position = [97 8 150 23];
            app.GenerateSummaryPlotsButton.Text = 'Generate Summary Plots';

            % Create FileOptionsPanel
            app.FileOptionsPanel = uipanel(app.SummaryPlotUIFigure);
            app.FileOptionsPanel.Title = 'File Options';
            app.FileOptionsPanel.Position = [10 293 328 113];

            % Create SelectOutputFolderButton
            app.SelectOutputFolderButton = uibutton(app.FileOptionsPanel, 'push');
            app.SelectOutputFolderButton.ButtonPushedFcn = createCallbackFcn(app, @SelectOutputFolderButtonPushed, true);
            app.SelectOutputFolderButton.WordWrap = 'on';
            app.SelectOutputFolderButton.Position = [12 43 88 41];
            app.SelectOutputFolderButton.Text = 'Select Output Folder';

            % Create OutputPathField
            app.OutputPathField = uieditfield(app.FileOptionsPanel, 'text');
            app.OutputPathField.Editable = 'off';
            app.OutputPathField.FontSize = 11;
            app.OutputPathField.Position = [111 54 202 22];

            % Create FigureFileNameEditFieldLabel
            app.FigureFileNameEditFieldLabel = uilabel(app.FileOptionsPanel);
            app.FigureFileNameEditFieldLabel.HorizontalAlignment = 'right';
            app.FigureFileNameEditFieldLabel.Position = [8 12 97 22];
            app.FigureFileNameEditFieldLabel.Text = 'Figure File Name';

            % Create FigureFileNameField
            app.FigureFileNameField = uieditfield(app.FileOptionsPanel, 'text');
            app.FigureFileNameField.ValueChangedFcn = createCallbackFcn(app, @FigureFileNameFieldValueChanged, true);
            app.FigureFileNameField.Position = [111 12 202 22];

            % Create BandLabelsPanel
            app.BandLabelsPanel = uipanel(app.SummaryPlotUIFigure);
            app.BandLabelsPanel.Title = 'Band Labels';
            app.BandLabelsPanel.Position = [10 38 328 244];

            % Create BandLabelsCommaSeparetedEditFieldLabel
            app.BandLabelsCommaSeparetedEditFieldLabel = uilabel(app.BandLabelsPanel);
            app.BandLabelsCommaSeparetedEditFieldLabel.HorizontalAlignment = 'center';
            app.BandLabelsCommaSeparetedEditFieldLabel.WordWrap = 'on';
            app.BandLabelsCommaSeparetedEditFieldLabel.Enable = 'off';
            app.BandLabelsCommaSeparetedEditFieldLabel.Position = [10 157 112 56];
            app.BandLabelsCommaSeparetedEditFieldLabel.Text = 'Band Labels (Comma Separeted)';

            % Create BandLabelsCommaSeparetedEditField
            app.BandLabelsCommaSeparetedEditField = uieditfield(app.BandLabelsPanel, 'text');
            app.BandLabelsCommaSeparetedEditField.Enable = 'off';
            app.BandLabelsCommaSeparetedEditField.Position = [138 174 173 22];

            % Create ControlLanesCheckBox
            app.ControlLanesCheckBox = uicheckbox(app.BandLabelsPanel);
            app.ControlLanesCheckBox.ValueChangedFcn = createCallbackFcn(app, @ControlLanesCheckBoxValueChanged, true);
            app.ControlLanesCheckBox.Enable = 'off';
            app.ControlLanesCheckBox.Text = 'Control Lanes';
            app.ControlLanesCheckBox.Position = [14 135 97 22];

            % Create LaneNumbersCommaSeparatedEditFieldLabel
            app.LaneNumbersCommaSeparatedEditFieldLabel = uilabel(app.BandLabelsPanel);
            app.LaneNumbersCommaSeparatedEditFieldLabel.HorizontalAlignment = 'center';
            app.LaneNumbersCommaSeparatedEditFieldLabel.WordWrap = 'on';
            app.LaneNumbersCommaSeparatedEditFieldLabel.Enable = 'off';
            app.LaneNumbersCommaSeparatedEditFieldLabel.Position = [12 75 110 56];
            app.LaneNumbersCommaSeparatedEditFieldLabel.Text = 'Lane Numbers (Comma Separated)';

            % Create LaneNumbersCommaSeparatedEditField
            app.LaneNumbersCommaSeparatedEditField = uieditfield(app.BandLabelsPanel, 'text');
            app.LaneNumbersCommaSeparatedEditField.Enable = 'off';
            app.LaneNumbersCommaSeparatedEditField.Position = [138 93 173 22];

            % Create ControlLabelsCommaSeparatedEditFieldLabel
            app.ControlLabelsCommaSeparatedEditFieldLabel = uilabel(app.SummaryPlotUIFigure);
            app.ControlLabelsCommaSeparatedEditFieldLabel.HorizontalAlignment = 'center';
            app.ControlLabelsCommaSeparatedEditFieldLabel.WordWrap = 'on';
            app.ControlLabelsCommaSeparatedEditFieldLabel.Enable = 'off';
            app.ControlLabelsCommaSeparatedEditFieldLabel.Position = [17 51 118 56];
            app.ControlLabelsCommaSeparatedEditFieldLabel.Text = 'Control Labels (Comma Separated)';

            % Create ControlLabelsCommaSeparatedEditField
            app.ControlLabelsCommaSeparatedEditField = uieditfield(app.SummaryPlotUIFigure, 'text');
            app.ControlLabelsCommaSeparatedEditField.Enable = 'off';
            app.ControlLabelsCommaSeparatedEditField.Position = [148 68 175 22];

            % Show the figure after all components are created
            app.SummaryPlotUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SummaryPlotWindow_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.SummaryPlotUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.SummaryPlotUIFigure)
        end
    end
end