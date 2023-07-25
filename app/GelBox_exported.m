classdef GelBox_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GelBoxUIFigure               matlab.ui.Figure
        FileMenu                     matlab.ui.container.Menu
        LoadImageMenu                matlab.ui.container.Menu
        InvertImageMenu              matlab.ui.container.Menu
        LoadAnalysisMenu             matlab.ui.container.Menu
        SaveAnalysisMenu             matlab.ui.container.Menu
        ExportResultsMenu            matlab.ui.container.Menu
        DataAnalysisMenu             matlab.ui.container.Menu
        GelImageFileInformationMenu  matlab.ui.container.Menu
        SelectedBoxInformationMenu   matlab.ui.container.Menu
        SummaryPlotMenu              matlab.ui.container.Menu
        FittingPanel                 matlab.ui.container.Panel
        FittingParametersButton      matlab.ui.control.Button
        NumberofBandsDropDown        matlab.ui.control.DropDown
        NumberofBandsDropDownLabel   matlab.ui.control.Label
        BandRelativeArea_3           matlab.ui.control.NumericEditField
        BandRelativeAreaLabel_3      matlab.ui.control.Label
        BandArea_3                   matlab.ui.control.NumericEditField
        BandAreaLabel_3              matlab.ui.control.Label
        BandRelativeArea_2           matlab.ui.control.NumericEditField
        BandRelativeAreaLabel_2      matlab.ui.control.Label
        BandArea_2                   matlab.ui.control.NumericEditField
        BandAreaLabel_2              matlab.ui.control.Label
        BandRelativeArea_1           matlab.ui.control.NumericEditField
        BandRelativeAreaLabel_1      matlab.ui.control.Label
        BandArea_1                   matlab.ui.control.NumericEditField
        BandAreaLabel_1              matlab.ui.control.Label
        rsquaredField                matlab.ui.control.NumericEditField
        RsquaredLabel                matlab.ui.control.Label
        DrawFittingCheckBox          matlab.ui.control.CheckBox
        BackgroundCorrectedOpticalDensityLabel  matlab.ui.control.Label
        RawOpticalDensityLabel       matlab.ui.control.Label
        raw_density_fit              matlab.ui.control.UIAxes
        background_corrected_raw_density_fit  matlab.ui.control.UIAxes
        GelImagePanel                matlab.ui.container.Panel
        BoxSelectionDropDown         matlab.ui.control.DropDown
        BoxSelectionDropDownLabel    matlab.ui.control.Label
        NewBoxButton                 matlab.ui.control.Button
        DeleteBoxButton              matlab.ui.control.Button
        gel_image_axes               matlab.ui.control.UIAxes
        OpticalDensitiesPanel        matlab.ui.container.Panel
        BackgroundAreaField          matlab.ui.control.NumericEditField
        BackgroundAreaLabel          matlab.ui.control.Label
        TotalAreaField               matlab.ui.control.NumericEditField
        TotalAreaEditFieldLabel      matlab.ui.control.Label
        BackgroundCorrectedOpticalDensityLabel_2  matlab.ui.control.Label
        RawOpticalDensityLabel_2     matlab.ui.control.Label
        BoxZoomLabel                 matlab.ui.control.Label
        box_inset                    matlab.ui.control.UIAxes
        background_corrected_raw_density  matlab.ui.control.UIAxes
        raw_density                  matlab.ui.control.UIAxes
    end


    properties (Access = public)
        gel_data % Description
        fitting_options
        new_box = 0
        mode_updated = 0
        parameters_updated = 0
        box_changed = 0
        loaded_analysis = 0
        par_est_na = 0
    end

    properties (Access = private)
        ImageFileTextDialog % Description
        SelectedBoxTextDialog % Description
        FittingOptions % Description
        SummaryPlot
    end

    methods (Access = public)

        function UpdateDisplay(app)
            %             temp_strings = get(gui.fitting_mode,'String');
            %             app.gel_data.fitting_mode = temp_strings{get(gui.fitting_mode,'Value')};
            if (isfield(app.gel_data,'box_handle'))
                n = numel(app.gel_data.box_handle);
                t = size(app.gel_data.par_update,2);
                if t ~= n
                    for l = t+1:n
                        app.gel_data.par_update(l) = 0;
                    end
                end
                % Get selected box in control
                control_strings = app.BoxSelectionDropDown.Value;
                selected_box = str2num(control_strings);

                for i=1:n
                    p(i,1:4) = app.gel_data.box_handle(i).Position;
                end
                w = p(selected_box,3);
                h = p(selected_box,4);
                if ((w~=app.gel_data.old_width)|(h~=app.gel_data.old_height))
                    for i=1:n
                        p(i,3) = w;
                        p(i,4) = h;
                        app.gel_data.box_handle(i).Position = p(i,:);
                    end
                end

                % Store data in case we need to save it
                for i=1:n
                    app.gel_data.box_position(i,:) = app.gel_data.box_handle(i).Position;
                end

                % Store data for display
                d=[];
                for i=1:n

                    d.box(i).fitting_mode = str2num(app.NumberofBandsDropDown.Value);
                    num_of_bands = d.box(i).fitting_mode;
                    % Extract position
                    d.box(i).position = app.gel_data.box_handle(i).Position;

                    % Label it
                    set(app.gel_data.box_label(i),'String',sprintf('%.0f',i));
                    set(app.gel_data.box_label(i), ...
                        'Position',[d.box(i).position(1)+d.box(i).position(3) ...
                        d.box(i).position(2)-50]);

                    % Calculate profile
                    d.box(i).inset = imcrop(app.gel_data.im_data, ...
                        d.box(i).position);

                    m = imcomplement(d.box(i).inset);

                    x = flipud(mean(m,2));
                    y = 1:size(m,1);

                    x_back = linspace(x(1),x(end),numel(y));
                    box_no = str2num(app.BoxSelectionDropDown.Value);
                    
                    if app.loaded_analysis || app.parameters_updated || app.gel_data.par_update(i)
                    elseif app.new_box || app.mode_updated || app.par_est_na
                        par_est = EstimateFittingParameters(app,y,x, ...
                            x_back,num_of_bands);
                        fnames = fieldnames(par_est);
                        for k = 1:numel(fnames)
                            app.gel_data.par_est(i).(fnames{k}) = ...
                                [];
                            app.gel_data.par_est(i).(fnames{k}) = ...
                                par_est.(fnames{k});
                        end
                    end

                    switch num_of_bands
                        case 1
                            app.BandRelativeArea_1.Enable = 0;
                            [x_bands,x_fit,r_squared] = ...
                                FitGaussian(app,y,x,x_back,i);
                        case 2
                            if ~app.BandRelativeArea_1.Enable
                                app.BandRelativeArea_1.Enable = 1;
                            end
                            [x_bands,x_fit,r_squared] = ...
                                Fit2Gaussian(app,y,x,x_back,i);
                        case 3
                            if ~app.BandRelativeArea_1.Enable
                                app.BandRelativeArea_1.Enable = 1;
                            end
                            [x_bands,x_fit,r_squared] = ...
                                Fit3Gaussian(app,y,x,x_back,i);
                    end

                    d.box(i).total_area = trapz(y,x);
                    d.box(i).background_area = trapz(y,x_back);

                    for j = 1 : num_of_bands
                        d.box(i).band_area(j) = trapz(y,x_bands(j,:));
                    end

                    % Store data for later

                    app.gel_data.box_data(i) = d.box(i);
                    app.gel_data.summary(i).x = x;
                    app.gel_data.summary(i).y = y;
                    app.gel_data.summary(i).x_fit = x_fit;
                    app.gel_data.summary(i).x_back = x_back;
                    

                    if num_of_bands == 2
                        [~,peak_band_1] = max(x_bands(1,:));
                        [~,peak_band_2] = max(x_bands(2,:));

                        if peak_band_2 > peak_band_1

                            app.gel_data.summary(i).bottom = d.box(i).band_area(1);
                            app.gel_data.summary(i).top = d.box(i).band_area(2);
                            app.gel_data.summary(i).band_1 = x_bands(1,:);
                            app.gel_data.summary(i).band_2 = x_bands(2,:);
                        else
                            app.gel_data.summary(i).bottom = d.box(i).band_area(2);
                            app.gel_data.summary(i).top = d.box(i).band_area(1);
                            app.gel_data.summary(i).band_1 = x_bands(2,:);
                            app.gel_data.summary(i).band_2 = x_bands(1,:);
                        end
                    elseif num_of_bands == 3
                        [~,peak_band_1] = max(x_bands(1,:));
                        [~,peak_band_2] = max(x_bands(2,:));
                        [~,peak_band_3] = max(x_bands(3,:));

                        peak_band = [peak_band_1 peak_band_2 peak_band_3];

                        %sort
                        [~,sort_ix] = sort(peak_band);
                        app.gel_data.summary(i).bottom = d.box(i).band_area( ...
                            sort_ix(1));
                        app.gel_data.summary(i).middle = d.box(i).band_area( ...
                            sort_ix(2));
                        app.gel_data.summary(i).top = d.box(i).band_area( ...
                            sort_ix(3));
                        app.gel_data.summary(i).band_1 = x_bands( ...
                            sort_ix(1),:);
                        app.gel_data.summary(i).band_2 = x_bands( ...
                            sort_ix(2),:);
                        app.gel_data.summary(i).band_3 = x_bands( ...
                            sort_ix(3),:);
                    else
                        app.gel_data.summary(i).band_1 = x_bands(1,:);
                        app.gel_data.summary(i).bottom = d.box(i).band_area;
                    end

                    app.gel_data.summary(i).inset = d.box(i).inset;
                    app.gel_data.summary(i).r_squared = r_squared;

                    % Display
                    if (i==selected_box)
                        app.TotalAreaField.Value = d.box(i).total_area;
                        app.BackgroundAreaField.Value = d.box(i).background_area;

                        app.rsquaredField.Value = r_squared;
                        center_image_with_preserved_aspect_ratio( ...
                            d.box(i).inset, ...
                            app.box_inset);
                        cla(app.raw_density)
                        plot(app.raw_density,x,y,"Color",'k',"LineWidth",2)
                        hold(app.raw_density,"on")
                        plot(app.raw_density,x_back,y,'-.m',"LineWidth",2)
                        x_pow = ceil(log10(max(x)));
                        x_tick_rounder = 10^(x_pow - 1);
                        x_t_end = ceil(max(x)/x_tick_rounder)*x_tick_rounder;
                        x_t_mid = round(x_t_end/2);
                        x_ticks = [0 x_t_mid x_t_end];
                        xticks(app.raw_density,x_ticks)
                        app.raw_density.XAxis.Exponent = 0;
                        xlim(app.raw_density,[0 x_t_end]);
                        ylim(app.raw_density,[1 max(y)]);
                        % legend(app.raw_density,'','Baseline', ...
                        %     'Location','best')


                        xticks(app.raw_density_fit,x_ticks);
                        app.raw_density_fit.XAxis.Exponent = 0;

                        xlim(app.raw_density_fit,[0 x_t_end]);
                        ylim(app.raw_density_fit,[1 max(y)]);
                        % legend(app.raw_density)


                        cla(app.background_corrected_raw_density)
                        plot(app.background_corrected_raw_density, ...
                            x-x_back',y,'-.k',"LineWidth",2)
                        hold(app.background_corrected_raw_density,"on")
                        plot(app.background_corrected_raw_density, ...
                            zeros(1,numel(y)),y,'-.m',"LineWidth",2)
                        % legend(app.background_corrected_raw_density,'','Baseline', ...
                        %     'Location','best')
                        ylim(app.background_corrected_raw_density, ...
                            [1 max(y)]);

                        x_pow = ceil(log10(max(x-x_back')));
                        x_tick_rounder = 10^(x_pow - 1);
                        x_t_end = ceil(max(x-x_back')/x_tick_rounder)*x_tick_rounder;
                        if min(x-x_back') < 0
                            x_tick_rounder = -10^(x_pow - 2)*0.25;
                        else
                            x_tick_rounder = 10^(x_pow - 2)*0.25;
                        end
                        x_t_beginning = ceil(min(x-x_back')/x_tick_rounder)*x_tick_rounder;
                        x_t_mid = round(x_t_end/2);
                        x_ticks = [x_t_beginning x_t_mid x_t_end];

                        app.background_corrected_raw_density.XAxis.Exponent = 0;
                        xticks(app.background_corrected_raw_density,x_ticks)
                        xlim(app.background_corrected_raw_density,[x_t_beginning x_t_end])
                        ylim(app.background_corrected_raw_density,[1 max(y)]);

                        xticks(app.background_corrected_raw_density_fit,x_ticks)
                        app.background_corrected_raw_density_fit.XAxis.Exponent = 0;
                        xlim(app.background_corrected_raw_density_fit,[x_t_beginning x_t_end])
                        ylim(app.background_corrected_raw_density_fit,[1 max(y)]);

                        switch num_of_bands
                            case 1
                                color = {'r'};
                                app.BandArea_1.Value = ...
                                    app.gel_data.summary(i).bottom;
                            case 2
                                color = {'r','b'};
                                total_band = app.gel_data.summary(i).bottom ...
                                    + app.gel_data.summary(i).top;
                                app.BandArea_1.Value = ...
                                    app.gel_data.summary(i).bottom;

                                app.BandRelativeArea_1.Value = ...
                                    app.gel_data.summary(i).bottom/total_band;

                                app.BandArea_2.Value = ...
                                    app.gel_data.summary(i).top;

                                app.BandRelativeArea_2.Value = ...
                                    app.gel_data.summary(i).top/total_band;
                            case 3
                                color = {'r','b','y'};
                                total_band = app.gel_data.summary(i).bottom ...
                                    + app.gel_data.summary(i).middle ...
                                    + app.gel_data.summary(i).top;

                                app.BandArea_1.Value = ...
                                    app.gel_data.summary(i).bottom;

                                app.BandRelativeArea_1.Value = ...
                                    app.gel_data.summary(i).bottom/total_band;

                                app.BandArea_2.Value = ...
                                    app.gel_data.summary(i).middle;

                                app.BandRelativeArea_2.Value = ...
                                    app.gel_data.summary(i).middle/total_band;

                                app.BandArea_3.Value = ...
                                    app.gel_data.summary(i).top;

                                app.BandRelativeArea_3.Value = ...
                                    app.gel_data.summary(i).top/total_band;
                        end


                        cla(app.raw_density_fit)
                        cla(app.background_corrected_raw_density_fit)



                        for j = 1 : num_of_bands
                            patch(app.raw_density_fit, ...
                                x_back+x_bands(j,:), ...
                                y,color{j},'FaceAlpha',0.65, ...
                                'EdgeColor',color{j},'EdgeAlpha',0.25, ...
                                'LineWidth',2)
                            patch(app.background_corrected_raw_density_fit, ...
                                x_bands(j,:), ...
                                y,color{j},'FaceAlpha',0.65, ...
                                'EdgeColor',color{j},'EdgeAlpha',0.25, ...
                                'LineWidth',2)
                        end

                        f_color = '#1aff00';

                        hold(app.background_corrected_raw_density_fit,"on")
                        plot(app.background_corrected_raw_density_fit, ...
                            x-x_back',y,'-.k',"LineWidth",2)
                        plot(app.background_corrected_raw_density_fit, ...
                            x_fit,y,':',"LineWidth",2,"Color",f_color)

                        hold(app.raw_density_fit,"on")

                        plot(app.raw_density_fit, ...
                            x,y,"Color",'k',"LineWidth",2)
                        plot(app.raw_density_fit, ...
                            x_back+x_fit,y,':',"LineWidth",2,"Color",f_color)

                        ylim(app.raw_density_fit, ...
                            [1 max(y)]);
                        ylim(app.background_corrected_raw_density_fit, ...
                            [1 max(y)]);
                    end
                end
                app.gel_data.d_box = d;
            end
        end

        function [y_bands, y_fit,r_squared] = FitGaussian(app,x,y,y_back,box_no)

            par_est = struct();
            box_pars = struct();
            box_pars = app.gel_data.par_est(box_no);
            names = fieldnames(box_pars);
            for m = 1 : numel(names)
                par_est.(names{m}) = box_pars.(names{m});
            end
            
            par = [par_est.peak_location(1) ...
                par_est.shape_parameter(1) ...
                par_est.amplitude(1) ...
                par_est.skew_parameter(1) ...
                ];
            
            target = y';

            target = target - y_back;

            j = 1;
            e = [];

            opts=optimset('fminsearch');
            opts.Display='off';
            opts.MaxIter=1000;
            opts.MaxFunEvals=10000;

            [p_result,fval,exitflag,output] = fminsearch(@profile_error_1gaussian, par, opts);

            r_squared = calculate_r_squared(y',y_fit+y_back);

            function trial_e = profile_error_1gaussian(par)

                [y_bands,y_fit] = calculate_profile(x,par);

                e(j) = 0;

                for i  = 1 : numel(target)
                    e(j) = e(j) + (y_fit(i) - target(i))^2;
                end

                trial_e = e(end);

                for i = 1:1
                    if any(y_bands(i,:)<0)
                        trial_e = trial_e + 10^12;
                    end
                end

                if any(par<0)
                    trial_e = trial_e + 10^12;
                end

                for i = 1:1
                    areas(i) = trapz(x,y_bands(i,:));
                end

                if any(areas<0)
                    trial_e = trial_e + 10^12;
                end

               positions = [par(1)];

                if any(positions>numel(x))
                    trial_e = trial_e + 10^12;
                end
                
                r_squared_iter = calculate_r_squared(target,y_fit);
                j = j + 1;
                if app.DrawFittingCheckBox.Value
                    cla(app.background_corrected_raw_density_fit)
                    plot(app.background_corrected_raw_density_fit,y_fit,x,'LineWidth',2)
                    hold(app.background_corrected_raw_density_fit, 'on')
                    plot(app.background_corrected_raw_density_fit,target,x,'LineWidth',2,'LineStyle','-.','Color','k')
                    app.rsquaredField.Value = r_squared_iter;
                    drawnow
                end
            end
            function [y_bands,y_fit] = calculate_profile(x,par)

                x1 = par(1);
                curve_shape1 = par(2);
                amp1 = par(3);
                skew1 = par(4);

                y_first = skewed_Gaussian(x,x1,curve_shape1,amp1,skew1);

                y_fit = y_first;
                y_bands(1,:) = y_first;

            end
            function y=skewed_Gaussian(x,x0,gamma,A,skew1)
                offset = zeros(1,length(x));
                offset((x-x0)>0) = skew1*(x((x-x0)>0)-x0);
                y=  A*exp(-gamma*(((x-x0)+offset).^2));
            end


        end

        function [y_bands, y_fit,r_squared] = Fit2Gaussian(app,x,y,...
                y_back,box_no)
            par_est = struct();
            box_pars = struct();
            box_pars = app.gel_data.par_est(box_no);
            names = fieldnames(box_pars);
            for m = 1 : numel(names)
                par_est.(names{m}) = box_pars.(names{m});
            end

            par = [par_est.peak_location(1) ...
                   par_est.shape_parameter(1)...
                   par_est.amplitude(1) ...
                   par_est.skew_parameter(1) ...
                   par_est.peak_location(2) ...
                   par_est.amplitude(2) ...
                  ];

            if numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1
                par(7) = par_est.shape_parameter(2);
                par(8) = par_est.skew_parameter(2);
            elseif numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) == 1
                par(7) = par_est.shape_parameter(2);
            elseif numel(unique(par_est.shape_parameter)) == 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1
                par(7) = par_est.skew_parameter(2);
            end


            j = 1;
            e = [];


            target = y';

            target = target - y_back;

            opts=optimset('fminsearch');
            opts.Display='off';
            opts.MaxIter=1000;
            opts.MaxFunEvals=10000;

            [p_result,fval,exitflag,output] = fminsearch(@profile_error_2gaussian, par, opts);

            p_result;
            r_squared = calculate_r_squared(y',y_fit+y_back);
            function trial_e = profile_error_2gaussian(par)

                [y_bands,y_fit] = calculate_2profile(x,par);

                e(j) = 0;

                for i  = 1 : numel(target)
                    e(j) = e(j) + (y_fit(i) - target(i))^2;
                end

                trial_e = e(end);

                for i = 1:2
                    if any(y_bands(i,:)<0)
                        trial_e = trial_e + 10^12;
                    end
                end

                if any(par<0)
                    trial_e = trial_e + 10^12;
                end

                for i = 1:2
                    areas(i) = trapz(x,y_bands(i,:));
                end

                if any(areas<0)
                    trial_e = trial_e + 10^12;
                end

               positions = [par(1) par(5)];

                if any(positions>numel(x))
                    trial_e = trial_e + 10^12;
                end

                r_squared_iter = calculate_r_squared(target,y_fit);
                j = j + 1;
                if app.DrawFittingCheckBox.Value
                    cla(app.background_corrected_raw_density_fit)
                    plot(app.background_corrected_raw_density_fit, ...
                        y_fit,x,'LineWidth',2)
                    hold(app.background_corrected_raw_density_fit, 'on')
                    plot(app.background_corrected_raw_density_fit, ...
                        target,x,'LineWidth',2,'LineStyle','-.','Color','k')
                    app.rsquaredField.Value = r_squared_iter;
                    drawnow
                end
            end
            function [y_bands,y_fit] = calculate_2profile(x,par)

                x1 = par(1);
                curve_shape1 = par(2);
                amp1 = par(3);
                skew1 = par(4);
                x2 = par(5);
                amp2 = par(6);

                if numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1
                    curve_shape2 = par(7);
                    skew2 = par(8);
                elseif numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) == 1
                    curve_shape2 = par(7);
                    skew2 = skew1;
                elseif numel(unique(par_est.shape_parameter)) == 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1
                    curve_shape2 = curve_shape1;
                    skew2 = par(7);
                else
                    curve_shape2 = curve_shape1;
                    skew2 = skew1;
                end

                y_first = skewed_Gaussian(x,x1,curve_shape1,amp1,skew1);
                y_second = skewed_Gaussian(x,x2,curve_shape2,amp2,skew2);

                y_fit = y_first + y_second;
                y_bands(1,:) = y_first;
                y_bands(2,:) = y_second;

            end
            function y=skewed_Gaussian(x,x0,gamma,A,skew1)
                offset = zeros(1,length(x));
                offset((x-x0)>0) = skew1*(x((x-x0)>0)-x0);
                y=  A*exp(-gamma*(((x-x0)+offset).^2));
            end
        end

        function [y_bands, y_fit,r_squared] = Fit3Gaussian(app,x,y,y_back,box_no)
            par_est = struct();
            box_pars = struct();
            box_pars = app.gel_data.par_est(box_no);
            names = fieldnames(box_pars);
            for m = 1 : numel(names)
                par_est.(names{m}) = box_pars.(names{m});
            end


            par = [par_est.peak_location(1) ...
                par_est.shape_parameter(1)...
                par_est.amplitude(1) ...
                par_est.skew_parameter(1) ...
                par_est.peak_location(2) ...
                par_est.amplitude(2) ...
                par_est.peak_location(3) ...
                par_est.amplitude(3) ...
                ];

            if numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1
                par(9) = par_est.shape_parameter(2);
                par(10) = par_est.skew_parameter(2);
                par(11) = par_est.shape_parameter(3);
                par(12) = par_est.skew_parameter(3);
            elseif numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) == 1
                par(9) = par_est.shape_parameter(2);
                par(10) = par_est.shape_parameter(3);
            elseif numel(unique(par_est.shape_parameter)) == 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1
                par(9) = par_est.skew_parameter(2);
                par(10) = par_est.skew_parameter(3);

            end

            
            j = 1;
            e = [];
            

            target = y';

            target = target - y_back;

            opts=optimset('fminsearch');
            opts.Display='off';
            opts.MaxIter=1000;
            opts.MaxFunEvals=10000;

            [p_result,fval,exitflag,output] = fminsearch(@profile_error_3gaussian, par, opts);

            r_squared = calculate_r_squared(y',y_fit+y_back);
            function trial_e = profile_error_3gaussian(par)
                [y_bands,y_fit] = calculate_3profile(x,par);
                
                e(j) = 0;
                for i  = 1 : numel(target)
                    e(j) = e(j) + (y_fit(i) - target(i))^2;
                end

                trial_e = e(end);

                for i = 1:3
                    if any(y_bands(i,:)<0)
                        trial_e = trial_e + 10^12;
                    end
                end

                if any(par<0)
                    trial_e = trial_e + 10^12;
                end

                for i = 1:3
                    areas(i) = trapz(x,y_bands(i,:));
                end

                if any(areas<0)
                    trial_e = trial_e + 10^12;
                end

                positions = [par(1) par(5) par(7)];

                if any(positions>numel(x))
                    trial_e = trial_e + 10^12;
                end

                r_squared_iter = calculate_r_squared(target,y_fit);
                j = j + 1;
                if app.DrawFittingCheckBox.Value
                    cla(app.background_corrected_raw_density_fit)
                    plot(app.background_corrected_raw_density_fit,y_fit,x,'LineWidth',2)
                    hold(app.background_corrected_raw_density_fit, 'on')
                    plot(app.background_corrected_raw_density_fit,target,x,'LineWidth',2,'LineStyle','-.','Color','k')
                    app.rsquaredField.Value = r_squared_iter;
                    drawnow
                end
            end
            function [y_bands,y_fit] = calculate_3profile(x,par)

                x1 = par(1);
                curve_shape1 = par(2);
                amp1 = par(3);
                skew1 = par(4);
                x2 = par(5);
                amp2 = par(6);
                x3 = par(7);
                amp3 = par(8);

                if numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1
                    curve_shape2 = par(9);
                    skew2 = par(10);
                    curve_shape3 = par(11);
                    skew3 = par(12);
                elseif numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) == 1
                    curve_shape2 = par(9);
                    skew2 = skew1;
                    curve_shape3 = par(10);
                    skew3 = skew1;
                elseif numel(unique(par_est.shape_parameter)) == 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1
                    curve_shape2 = curve_shape1;
                    skew2 = par(9);
                    curve_shape3 = curve_shape1;
                    skew3 = par(10);
                else
                    curve_shape2 = curve_shape1;
                    skew2 = skew1;
                    curve_shape3 = curve_shape1;
                    skew3 = skew1;
                end

                y_first = skewed_Gaussian(x,x1,curve_shape1,amp1,skew1);
                y_second = skewed_Gaussian(x,x2,curve_shape2,amp2,skew2);
                y_third = skewed_Gaussian(x,x3,curve_shape3,amp3,skew3);


                y_fit = y_first + y_second + y_third;
                y_bands(1,:) = y_first;
                y_bands(2,:) = y_second;
                y_bands(3,:) = y_third;

            end
            function y=skewed_Gaussian(x,x0,gamma,A,skew1)
                offset = zeros(1,length(x));
                offset((x-x0)>0) = skew1*(x((x-x0)>0)-x0);
                y=  A*exp(-gamma*(((x-x0)+offset).^2));
            end
        end

        function UpdateFittingOptions(app)
            app.parameters_updated = 1;
            UpdateDisplay(app)
            app.parameters_updated = 0;
        end
        
        function par_est = EstimateFittingParameters(app,x,y,y_back, ...
                no_of_bands)

            switch no_of_bands
                case 1 

                    peaks=find_peaks('x',x, ...
                        'y',y, ...
                        'min_rel_delta_y',0.05, ...
                        'min_x_index_spacing',2);

                     if numel(peaks.max_indices) == no_of_bands
                         first_curve_x_estimate=peaks.max_indices(1);
                     else
                         first_curve_x_estimate = 0.5*length(x);
                     end

                     target = y';

                     target = target - y_back;
                     [max_value,~]=max(target);

                     half_distance=(0.1*length(x));
                     alfa_estimate = -log(0.5)/(half_distance^2);

                     first_curve_shape_estimate = alfa_estimate;
                     first_curve_amp_estimate = max_value;
                     first_curve_skew_estimate = 1;


                     par_est.band_no  = 1;
                     par_est.peak_location = first_curve_x_estimate;
                     par_est.shape_parameter = first_curve_shape_estimate;
                     par_est.amplitude = first_curve_amp_estimate;
                     par_est.skew_parameter = first_curve_skew_estimate;

                case 2

                    peaks=find_peaks('x',x, ...
                        'y',y, ...
                        'min_rel_delta_y',0.05, ...
                        'min_x_index_spacing',2);

                    if numel(peaks.max_indices) == no_of_bands
                        first_curve_x_estimate=peaks.max_indices(1);
                        second_curve_x_estimate=peaks.max_indices(2);
                    else
                        first_curve_x_estimate = 0.3*length(x);
                        second_curve_x_estimate = 0.6*length(x);
                    end

                    target = y';

                    target = target - y_back;
                    [max_value,~]=max(target);

                    half_distance=(0.1*length(x));
                    alfa_estimate = -log(0.5)/(half_distance^2);

                    first_curve_shape_estimate = alfa_estimate;
                    first_curve_amp_estimate = max_value;
                    first_curve_skew_estimate = 1;

                    second_curve_amp_estimate = first_curve_amp_estimate;
                    second_curve_shape_estimate = alfa_estimate;
                    second_curve_skew_estimate = 1;

                     par_est.band_no  = [1;2];
                     par_est.peak_location = ...
                     [first_curve_x_estimate;second_curve_x_estimate];
                     par_est.shape_parameter = ...
                     [first_curve_shape_estimate;second_curve_shape_estimate];
                     par_est.amplitude = ...
                         [first_curve_amp_estimate;second_curve_amp_estimate];
                     par_est.skew_parameter = ...
                     [first_curve_skew_estimate;second_curve_skew_estimate];

                case 3
                    peaks=find_peaks('x',x, ...
                        'y',y, ...
                        'min_rel_delta_y',0.05, ...
                        'min_x_index_spacing',2);
                    if numel(peaks.max_indices) == no_of_bands
                        first_curve_x_estimate=peaks.max_indices(1);
                        second_curve_x_estimate=peaks.max_indices(2);
                        third_curve_x_estimate=peaks.max_indices(3);
                    else
                        first_curve_x_estimate = 0.2*length(x);
                        second_curve_x_estimate = 0.3*length(x);
                        third_curve_x_estimate = 0.6*length(x);
                    end
                    target = y';

                    target = target - y_back;
                    [max_value,~]=max(target);

                    half_distance=(0.1*length(x));
                    alfa_estimate = -log(0.5)/(half_distance^2);

                    first_curve_shape_estimate = alfa_estimate;
                    first_curve_amp_estimate = max_value;
                    first_curve_skew_estimate = 1;

                    second_curve_amp_estimate = first_curve_amp_estimate;
                    second_curve_shape_estimate = alfa_estimate;
                    second_curve_skew_estimate = 1;

                    third_curve_amp_estimate = first_curve_amp_estimate;
                    third_curve_shape_estimate = alfa_estimate;
                    third_curve_skew_estimate = 1;

                    par_est.band_no  = [1;2;3];
                    par_est.peak_location = ...
                        [first_curve_x_estimate;...
                        second_curve_x_estimate;...
                        third_curve_x_estimate];
                    par_est.shape_parameter = ...
                        [first_curve_shape_estimate;...
                        second_curve_shape_estimate;...
                        third_curve_shape_estimate];
                    par_est.amplitude = ...
                        [first_curve_amp_estimate;...
                        second_curve_amp_estimate;...
                        third_curve_amp_estimate];
                    par_est.skew_parameter = ...
                        [first_curve_skew_estimate;...
                        second_curve_skew_estimate;...
                        third_curve_skew_estimate];
            end


            
        end
    end

    methods (Access = private)

        function ResetDisplay(app)

            % Clear the figure displays
            cla(app.box_inset)
            cla(app.raw_density)
            cla(app.background_corrected_raw_density)
            cla(app.raw_density_fit)
            cla(app.background_corrected_raw_density_fit)

            % Reset optical density fields
            app.TotalAreaField.Value = 0;
            app.BackgroundAreaField.Value = 0;

            % Reset fitting fields
            app.DrawFittingCheckBox.Value = 0;
            app.rsquaredField.Value = 0;
            app.BandArea_1.Value = 0;
            app.BandArea_2.Value = 0;
            app.BandArea_3.Value = 0;
            app.BandRelativeArea_1.Value = 0;
            app.BandRelativeArea_2.Value = 0;
            app.BandRelativeArea_3.Value = 0;

            % Reset number of bands
            app.NumberofBandsDropDown.Value = '1';

            % Reset box controls
            app.DeleteBoxButton.Enable = 0;
            control_strings = {'1'};
            app.BoxSelectionDropDown.Items = control_strings;
            app.BoxSelectionDropDown.Value = control_strings{1};

        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)

            addpath(genpath('../code/utilities'));
            movegui(app.GelBoxUIFigure,'center')
            
            
            
            colormap(app.GelBoxUIFigure, 'gray');
            app.fitting_options.shared_shape = true(1,1);
            app.fitting_options.shared_skewness = true(1,1);
            disableDefaultInteractivity(app.gel_image_axes)

        end

        % Menu selected function: LoadImageMenu
        function LoadImageButtonPushed(app, event)
            [file_string,path_string]=uigetfile2( ...
                {'*.png','PNG';'*.tif','TIF'}, ...
                'Select image file');
            if (path_string~=0)

                ResetDisplay(app)
                app.DataAnalysisMenu.Enable = 1;
                app.gel_data = [];
                app.gel_data.par_update = 0;

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

            end

        end

        % Menu selected function: InvertImageMenu
        function InvertImageButtonPushed(app, event)
            app.gel_data.im_data = imcomplement(app.gel_data.im_data);
            center_image_with_preserved_aspect_ratio( ...
                app.gel_data.im_data, ...
                app.gel_image_axes);
            app.gel_data.invert_status = 1;
        end

        % Menu selected function: LoadAnalysisMenu
        function LoadAnalysisButtonPushed(app, event)

            [file_string,path_string] = uigetfile2( ...
                {'*.gdf','Gel data file'},'Select file to load prior analysis');

            if (path_string~=0)
            % Delete any old boxes
            if (isfield(app.gel_data,'box_handle'))
                n = numel(app.gel_data.box_handle);
                for i=1:n
                    delete(app.gel_data.box_handle(i));
                end
            end
                app.DeleteBoxButton.Enable = 1;
                app.DataAnalysisMenu.Enable = 1;
                
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

                n=size(save_data.box_position,1);
                names = {'band_no','peak_location','shape_parameter','amplitude','skew_parameter'};
                for j = 1:n
                    for i = 1:numel(names)
                        if isfield(save_data,'par_est')
                        app.gel_data.par_est(j).(names{i}) = save_data.par_est(j).(names{i});
                        app.loaded_analysis = 1;
                        else
                        app.gel_data.par_est(j).(names{i}) = [];
                        app.par_est_na = 1;
                        end
                    end
                end
                control_strings = [];
                app.gel_data.par_update = zeros(1,n);
                
                if isfield(save_data,'par_est') && app.gel_data.par_est(1).band_no(end) ~= 1
                    val = app.gel_data.par_est(1).band_no(end);
                    app.NumberofBandsDropDown.Value = num2str(val);
                    switch val
                        case 1
                            try
                                app.BandRelativeAreaLabel_1.Value = 0;
                                app.BandRelativeAreaLabel_1.Enable = 0;
                                app.BandArea_2.Value = 0;
                                app.BandArea_2.Enable = 0;
                                app.BandRelativeArea_2.Value = 0;
                                app.BandRelativeArea_2.Enable = 0;
                                app.BandAreaLabel_2.Enable = 0;
                                app.BandRelativeAreaLabel_2.Enable = 0;
                                app.BandArea_3.Value = 0;
                                app.BandArea_3.Enable = 0;
                                app.BandRelativeArea_3.Value = 0;
                                app.BandRelativeArea_3.Enable = 0;
                                app.BandAreaLabel_3.Enable = 0;
                                app.BandRelativeAreaLabel_3.Enable = 0;
                            end
                        case 2
                            app.BandRelativeArea_1.Enable = 1;
                            app.BandRelativeAreaLabel_1.Enable = 1;
                            app.BandArea_2.Enable = 1;
                            app.BandRelativeArea_2.Enable = 1;
                            app.BandAreaLabel_2.Enable = 1;
                            app.BandRelativeAreaLabel_2.Enable = 1;
                            try
                                app.BandArea_3.Value = 0;
                                app.BandArea_3.Enable = 0;
                                app.BandRelativeArea_3.Enable = 0;
                                app.BandRelativeArea_3.Value = 0;
                                app.BandAreaLabel_3.Enable = 0;
                                app.BandRelativeAreaLabel_3.Enable = 0;
                            end

                        case 3
                            app.BandRelativeArea_1.Enable = 1;
                            app.BandRelativeAreaLabel_1.Enable = 1;
                            app.BandArea_2.Enable = 1;
                            app.BandRelativeArea_2.Enable = 1;
                            app.BandAreaLabel_2.Enable = 1;
                            app.BandRelativeAreaLabel_2.Enable = 1;
                            app.BandArea_3.Enable = 1;
                            app.BandRelativeArea_3.Enable = 1;
                            app.BandAreaLabel_3.Enable = 1;
                            app.BandRelativeAreaLabel_3.Enable = 1;
                    end
                end

                for i=1:n
                    app.gel_data.box_handle(i) = images.roi.Rectangle(app.gel_image_axes, ...
                        'Position',save_data.box_position(i,:));
                    control_strings{i} = sprintf('%.0f',i);
                end

                app.BoxSelectionDropDown.Items = control_strings;
                app.BoxSelectionDropDown.Value = control_strings{1};

                for i=1:n
                    app.gel_data.box_handle(i).FaceAlpha = 0;
                    if (i~=1)
                        app.gel_data.box_handle(i).Color = [1 0 0];
                        app.gel_data.box_handle(i).InteractionsAllowed = 'none';
                    else
                        app.gel_data.box_handle(i).Color = [0 1 0];
                        app.gel_data.box_handle(i).InteractionsAllowed = 'all';

                    end

                    p = app.gel_data.box_handle(i).Position;
                    app.gel_data.box_label(i) = text(p(1)+p(3),p(2)-50,sprintf('%.0f',i), ...
                        'Parent',app.gel_image_axes);

                    app.gel_data.old_width = p(3);
                    app.gel_data.old_height = p(4);

                    i=i;
                    addlistener(app.gel_data.box_handle(i),"MovingROI",@(src,evt) new_box_position2(evt));
                end

                % Need this to make labels
                drawnow;
            end

            UpdateDisplay(app)
            app.loaded_analysis = 0;
            app.par_est_na = 0;

            % Nested function
            function new_box_position2(evt);
                if (isfield(app.gel_data,'box_position'))
                    box_position = app.gel_data.box_position;
                    [r,c]=size(box_position);
                    if (r>=n)&(~isequal(box_position(n,:),evt.CurrentPosition))
                        app.new_box = n;
                        UpdateDisplay(app)
                    end
                else
                    app.new_box = n;
                    UpdateDisplay(app)
                end
            end
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
                    'Position',p + [120,0,0,0]);
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

            app.DeleteBoxButton.Enable = 1;
            app.SelectedBoxInformationMenu.Enable = 1;
            app.new_box = n;
            UpdateDisplay(app);

            function new_box_position(evt);
                if (isfield(app.gel_data,'box_position'))
                    box_position = app.gel_data.box_position;
                    [r,c]=size(box_position);
                    if (r>=n)&(~isequal(box_position(n,:),evt.CurrentPosition))
                        app.new_box = n;
                        UpdateDisplay(app);
                    end
                else
                    app.new_box = n;
                    UpdateDisplay(app);
                end
            end

        end

        % Value changed function: NumberofBandsDropDown
        function NumberofBandsDropDownValueChanged(app, event)
            value = app.NumberofBandsDropDown.Value;
            app.mode_updated = 1;

            switch value
                case '1'
                    try
                        app.BandRelativeAreaLabel_1.Enable = 0;
                        app.BandArea_2.Value = 0;
                        app.BandArea_2.Enable = 0;
                        app.BandRelativeArea_2.Value = 0;
                        app.BandRelativeArea_2.Enable = 0;
                        app.BandAreaLabel_2.Enable = 0;
                        app.BandRelativeAreaLabel_2.Enable = 0;
                        app.BandArea_3.Value = 0;
                        app.BandArea_3.Enable = 0;
                        app.BandRelativeArea_3.Value = 0;
                        app.BandRelativeArea_3.Enable = 0;
                        app.BandAreaLabel_3.Enable = 0;
                        app.BandRelativeAreaLabel_3.Enable = 0;
                    end
                case '2'
                    app.BandRelativeArea_1.Enable = 1;
                    app.BandRelativeAreaLabel_1.Enable = 1;
                    app.BandArea_2.Enable = 1;
                    app.BandRelativeArea_2.Enable = 1;
                    app.BandAreaLabel_2.Enable = 1;
                    app.BandRelativeAreaLabel_2.Enable = 1;
                    try
                        app.BandArea_3.Value = 0;
                        app.BandArea_3.Enable = 0;
                        app.BandRelativeArea_3.Enable = 0;
                        app.BandRelativeArea_3.Value = 0;
                        app.BandAreaLabel_3.Enable = 0;
                        app.BandRelativeAreaLabel_3.Enable = 0;
                    end

                case '3'
                    app.BandRelativeArea_1.Enable = 1;
                    app.BandRelativeAreaLabel_1.Enable = 1;
                    app.BandArea_2.Enable = 1;
                    app.BandRelativeArea_2.Enable = 1;
                    app.BandAreaLabel_2.Enable = 1;
                    app.BandRelativeAreaLabel_2.Enable = 1;
                    app.BandArea_3.Enable = 1;
                    app.BandRelativeArea_3.Enable = 1;
                    app.BandAreaLabel_3.Enable = 1;
                    app.BandRelativeAreaLabel_3.Enable = 1;
            end
            UpdateDisplay(app);
            app.mode_updated = 0;

        end

        % Menu selected function: SaveAnalysisMenu
        function SaveAnalysisButtonPushed(app, event)
            save_data.image_file_string = app.gel_data.image_file_string;
            try
                save_data.box_position = app.gel_data.box_position;
            catch
                warndlg('The analysis boxes are not available.')
                return
            end
            save_data.im_data = app.gel_data.im_data;
            save_data.imfinfo = app.gel_data.imfinfo;
            
            number_of_boxes = size(save_data.box_position,1);
            names = {'band_no','peak_location','shape_parameter','amplitude','skew_parameter'};
            for j = 1:number_of_boxes
                for i = 1:numel(names)
                    save_data.par_est(j).(names{i}) = app.gel_data.par_est(j).(names{i});
                end
            end

            [file_string,path_string] = uiputfile2( ...
                {'*.gdf','Gel data file'},'Select file to save analysis');

            if (path_string~=0)
                save(fullfile(path_string,file_string),'save_data');

                msgbox(sprintf('Current analysis saved to %s',file_string), ...
                    'Analysis saved');
            end
        end

        % Value changed function: BoxSelectionDropDown
        function BoxSelectionDropDownValueChanged(app, event)
            selected_box = str2num(app.BoxSelectionDropDown.Value);
            control_strings = app.BoxSelectionDropDown.Items;

            n = numel(control_strings);
            for i=1:n
                if (i~=selected_box)
                    app.gel_data.box_handle(i).Color = [1 0 0];
                    app.gel_data.box_handle(i).InteractionsAllowed = 'none';
                else
                    app.gel_data.box_handle(i).Color = [0 1 0];
                    app.gel_data.box_handle(i).InteractionsAllowed = 'all';
                end
            end
%             app.new_box = 0;
            UpdateDisplay(app)


        end

        % Button pushed function: DeleteBoxButton
        function DeleteBoxButtonPushed(app, event)
            selected_box = str2num(app.BoxSelectionDropDown.Value);
            control_strings = app.BoxSelectionDropDown.Items;

            delete(app.gel_data.box_handle(selected_box))
            delete(app.gel_data.box_label(selected_box))
            app.gel_data.box_handle(selected_box) = [];
            app.gel_data.box_label(selected_box) = [];

            n = numel(control_strings) - 1;
            control_strings = {};
            for i=1:n
                control_strings{i} = sprintf('%.0f',i);
            end

            app.BoxSelectionDropDown.Items = control_strings;
            if selected_box == 1
                app.BoxSelectionDropDown.Value = num2str(selected_box);
            else
                app.BoxSelectionDropDown.Value = num2str(selected_box-1);
            end

            new_selected_box = app.BoxSelectionDropDown.Value;

            for i=1:n
                delete(app.gel_data.box_label(i));
                if (i~=new_selected_box)
                    app.gel_data.box_handle(i).Color = [1 0 0];
                    app.gel_data.box_handle(i).InteractionsAllowed = 'none';
                else
                    app.gel_data.box_handle(i).Color = [0 1 0];
                    app.gel_data.box_handle(i).InteractionsAllowed = 'all';
                end
            end

            app.gel_data.box_label = [];
            for i = 1:n
                p = app.gel_data.box_handle(i).Position;
                app.gel_data.box_label(i) = text(p(1)+p(3),p(2)-50,sprintf('%.0f',i), ...
                    'Parent',app.gel_image_axes);

                app.gel_data.old_width = p(3);
                app.gel_data.old_height = p(4);
            end
        end

        % Menu selected function: ExportResultsMenu
        function OutputButtonPushed(app, event)
            d = [];
            d.image_file{1} = app.gel_data.image_file_string;
            n = numel(app.gel_data.box_handle);
            for i = 2 : n
            d.image_file{i,1} = '';
            end
            for i=1:n
                d.box(i) = i;
                d.total_area(i) = app.gel_data.box_data(i).total_area;
                d.background_area(i) = app.gel_data.box_data(i).background_area;
                num_of_bands = numel(app.gel_data.box_data(i).band_area);
                switch num_of_bands
                    case 1
                        d.band_area_bottom(i) = app.gel_data.summary(i).bottom;
                    case 2
                        d.band_area_bottom(i) = app.gel_data.summary(i).bottom;
                        d.band_area_top(i) = app.gel_data.summary(i).top;
                    case 3
                        d.band_area_bottom(i) = app.gel_data.summary(i).bottom;
                        d.band_area_middle(i) = app.gel_data.summary(i).middle;
                        d.band_area_top(i) = app.gel_data.summary(i).top;
                end

                d.band_left(i) = app.gel_data.box_data(i).position(1);
                d.band_top(i) = app.gel_data.box_data(i).position(2);
                d.band_width(i) = app.gel_data.box_data(i).position(3);
                d.band_height(i) = app.gel_data.box_data(i).position(4);
                d.fitting_mode{i} = app.gel_data.box_data(i).fitting_mode;
                d.num_of_bands(i) = num_of_bands;
                d.r_squared(i) = app.gel_data.summary(i).r_squared;
            end

            [file_string,path_string] = uiputfile2( ...
                {'*.xlsx','Excel file'},'Select file for results');
            
            d_out = d;
            names = fieldnames(d_out);
            
            for j = 1 : length(d_out)
                for i = 1 : numel(names)

                    [row,col] = size(d_out(j).(names{i}));

                    if row == 1 && col ~= 1
                        d_out(j).(names{i}) = (d_out(j).(names{i}))';
                    end
                end
                
            end
            
            output_file = fullfile(path_string,file_string);
            
            try
                delete(output_file);
            end
            
            writetable(struct2table(d_out),output_file,'Sheet','Summary')
             
            summary = app.gel_data.summary;

            for i = 1:n
                i_name = sprintf('box_%i',i);
            end

            switch num_of_bands
                case 1
                    summary = rmfield(summary,'bottom');
                    try
                    summary = rmfield(summary,'middle');
                    summary = rmfield(summary,'top');
                    end
                case 2
                    summary = rmfield(summary,'bottom');
                    summary = rmfield(summary,'top');
                    try
                    summary = rmfield(summary,'middle');
                    end

                case 3
                    summary = rmfield(summary,'bottom');
                    summary = rmfield(summary,'middle');
                    summary = rmfield(summary,'top');

            end
            summary = rmfield(summary,'r_squared');
            summary = rmfield(summary,'inset');

            names = fieldnames(summary);
            for j = 1 : length(summary)
                for i = 1 : numel(names)

                    [row,col] = size(summary(j).(names{i}));

                    if row == 1 && col ~= 1
                        summary(j).(names{i}) = (summary(j).(names{i}))';
                    end
                end
                sheet = sprintf('box_%i',j);
                writetable(struct2table(summary(j)),output_file,'Sheet',sheet)
            end


        end

        % Menu selected function: GelImageFileInformationMenu
        function GelImageFileInformationMenuSelected(app, event)
            app.ImageFileTextDialog = ImFInfoDialog(app);
        end

        % Close request function: GelBoxUIFigure
        function GelBoxUIFigureCloseRequest(app, event)
            delete(app.ImageFileTextDialog)
            delete(app.SelectedBoxTextDialog)
            delete(app.FittingOptions)
            delete(app)

        end

        % Menu selected function: SelectedBoxInformationMenu
        function SelectedBoxInformationMenuSelected(app, event)
            app.SelectedBoxTextDialog = SelectedBoxDialog(app);
        end

        % Button pushed function: FittingParametersButton
        function FittingParametersButtonPushed(app, event)
            app.FittingOptions = FittingOptionsDialog(app);

        end

        % Menu selected function: SummaryPlotMenu
        function AnalysisSummaryPlotMenuSelected(app, event)
            app.SummaryPlot = SummaryPlotWindow(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create GelBoxUIFigure and hide until all components are created
            app.GelBoxUIFigure = uifigure('Visible', 'off');
            app.GelBoxUIFigure.Colormap = [0.2431 0.149 0.6588;0.2431 0.1529 0.6745;0.2471 0.1569 0.6863;0.2471 0.1608 0.698;0.251 0.1647 0.7059;0.251 0.1686 0.7176;0.2549 0.1725 0.7294;0.2549 0.1765 0.7412;0.2588 0.1804 0.749;0.2588 0.1804 0.7608;0.2627 0.1882 0.7725;0.2588 0.1882 0.7804;0.2627 0.1961 0.7922;0.2667 0.2 0.8039;0.2667 0.2039 0.8157;0.2706 0.2078 0.8235;0.2706 0.2157 0.8353;0.2706 0.2196 0.8431;0.2745 0.2235 0.851;0.2745 0.2275 0.8627;0.2745 0.2314 0.8706;0.2745 0.2392 0.8784;0.2784 0.2431 0.8824;0.2784 0.2471 0.8902;0.2784 0.2549 0.898;0.2784 0.2588 0.902;0.2784 0.2667 0.9098;0.2784 0.2706 0.9137;0.2784 0.2745 0.9216;0.2824 0.2824 0.9255;0.2824 0.2863 0.9294;0.2824 0.2941 0.9333;0.2824 0.298 0.9412;0.2824 0.3059 0.9451;0.2824 0.3098 0.949;0.2824 0.3137 0.9529;0.2824 0.3216 0.9569;0.2824 0.3255 0.9608;0.2824 0.3294 0.9647;0.2784 0.3373 0.9686;0.2784 0.3412 0.9686;0.2784 0.349 0.9725;0.2784 0.3529 0.9765;0.2784 0.3569 0.9804;0.2784 0.3647 0.9804;0.2745 0.3686 0.9843;0.2745 0.3765 0.9843;0.2745 0.3804 0.9882;0.2706 0.3843 0.9882;0.2706 0.3922 0.9922;0.2667 0.3961 0.9922;0.2627 0.4039 0.9922;0.2627 0.4078 0.9961;0.2588 0.4157 0.9961;0.2549 0.4196 0.9961;0.251 0.4275 0.9961;0.2471 0.4314 1;0.2431 0.4392 1;0.2353 0.4431 1;0.2314 0.451 1;0.2235 0.4549 1;0.2196 0.4627 0.9961;0.2118 0.4667 0.9961;0.2078 0.4745 0.9922;0.2 0.4784 0.9922;0.1961 0.4863 0.9882;0.1922 0.4902 0.9882;0.1882 0.498 0.9843;0.1843 0.502 0.9804;0.1843 0.5098 0.9804;0.1804 0.5137 0.9765;0.1804 0.5176 0.9725;0.1804 0.5255 0.9725;0.1804 0.5294 0.9686;0.1765 0.5333 0.9647;0.1765 0.5412 0.9608;0.1765 0.5451 0.9569;0.1765 0.549 0.9529;0.1765 0.5569 0.949;0.1725 0.5608 0.9451;0.1725 0.5647 0.9412;0.1686 0.5686 0.9373;0.1647 0.5765 0.9333;0.1608 0.5804 0.9294;0.1569 0.5843 0.9255;0.1529 0.5922 0.9216;0.1529 0.5961 0.9176;0.149 0.6 0.9137;0.149 0.6039 0.9098;0.1451 0.6078 0.9098;0.1451 0.6118 0.9059;0.1412 0.6196 0.902;0.1412 0.6235 0.898;0.1373 0.6275 0.898;0.1373 0.6314 0.8941;0.1333 0.6353 0.8941;0.1294 0.6392 0.8902;0.1255 0.6471 0.8902;0.1216 0.651 0.8863;0.1176 0.6549 0.8824;0.1137 0.6588 0.8824;0.1137 0.6627 0.8784;0.1098 0.6667 0.8745;0.1059 0.6706 0.8706;0.102 0.6745 0.8667;0.098 0.6784 0.8627;0.0902 0.6824 0.8549;0.0863 0.6863 0.851;0.0784 0.6902 0.8471;0.0706 0.6941 0.8392;0.0627 0.698 0.8353;0.0549 0.702 0.8314;0.0431 0.702 0.8235;0.0314 0.7059 0.8196;0.0235 0.7098 0.8118;0.0157 0.7137 0.8078;0.0078 0.7176 0.8;0.0039 0.7176 0.7922;0 0.7216 0.7882;0 0.7255 0.7804;0 0.7294 0.7765;0.0039 0.7294 0.7686;0.0078 0.7333 0.7608;0.0157 0.7333 0.7569;0.0235 0.7373 0.749;0.0353 0.7412 0.7412;0.051 0.7412 0.7373;0.0627 0.7451 0.7294;0.0784 0.7451 0.7216;0.0902 0.749 0.7137;0.102 0.7529 0.7098;0.1137 0.7529 0.702;0.1255 0.7569 0.6941;0.1373 0.7569 0.6863;0.1451 0.7608 0.6824;0.1529 0.7608 0.6745;0.1608 0.7647 0.6667;0.1686 0.7647 0.6588;0.1725 0.7686 0.651;0.1804 0.7686 0.6471;0.1843 0.7725 0.6392;0.1922 0.7725 0.6314;0.1961 0.7765 0.6235;0.2 0.7804 0.6157;0.2078 0.7804 0.6078;0.2118 0.7843 0.6;0.2196 0.7843 0.5882;0.2235 0.7882 0.5804;0.2314 0.7882 0.5725;0.2392 0.7922 0.5647;0.251 0.7922 0.5529;0.2588 0.7922 0.5451;0.2706 0.7961 0.5373;0.2824 0.7961 0.5255;0.2941 0.7961 0.5176;0.3059 0.8 0.5059;0.3176 0.8 0.498;0.3294 0.8 0.4863;0.3412 0.8 0.4784;0.3529 0.8 0.4667;0.3686 0.8039 0.4549;0.3804 0.8039 0.4471;0.3922 0.8039 0.4353;0.4039 0.8039 0.4235;0.4196 0.8039 0.4118;0.4314 0.8039 0.4;0.4471 0.8039 0.3922;0.4627 0.8 0.3804;0.4745 0.8 0.3686;0.4902 0.8 0.3569;0.5059 0.8 0.349;0.5176 0.8 0.3373;0.5333 0.7961 0.3255;0.5451 0.7961 0.3176;0.5608 0.7961 0.3059;0.5765 0.7922 0.2941;0.5882 0.7922 0.2824;0.6039 0.7882 0.2745;0.6157 0.7882 0.2627;0.6314 0.7843 0.251;0.6431 0.7843 0.2431;0.6549 0.7804 0.2314;0.6706 0.7804 0.2235;0.6824 0.7765 0.2157;0.698 0.7765 0.2078;0.7098 0.7725 0.2;0.7216 0.7686 0.1922;0.7333 0.7686 0.1843;0.7451 0.7647 0.1765;0.7608 0.7647 0.1725;0.7725 0.7608 0.1647;0.7843 0.7569 0.1608;0.7961 0.7569 0.1569;0.8078 0.7529 0.1529;0.8157 0.749 0.1529;0.8275 0.749 0.1529;0.8392 0.7451 0.1529;0.851 0.7451 0.1569;0.8588 0.7412 0.1569;0.8706 0.7373 0.1608;0.8824 0.7373 0.1647;0.8902 0.7373 0.1686;0.902 0.7333 0.1765;0.9098 0.7333 0.1804;0.9176 0.7294 0.1882;0.9255 0.7294 0.1961;0.9373 0.7294 0.2078;0.9451 0.7294 0.2157;0.9529 0.7294 0.2235;0.9608 0.7294 0.2314;0.9686 0.7294 0.2392;0.9765 0.7294 0.2431;0.9843 0.7333 0.2431;0.9882 0.7373 0.2431;0.9961 0.7412 0.2392;0.9961 0.7451 0.2353;0.9961 0.7529 0.2314;0.9961 0.7569 0.2275;0.9961 0.7608 0.2235;0.9961 0.7686 0.2196;0.9961 0.7725 0.2157;0.9961 0.7804 0.2078;0.9961 0.7843 0.2039;0.9961 0.7922 0.2;0.9922 0.7961 0.1961;0.9922 0.8039 0.1922;0.9922 0.8078 0.1922;0.9882 0.8157 0.1882;0.9843 0.8235 0.1843;0.9843 0.8275 0.1804;0.9804 0.8353 0.1804;0.9765 0.8392 0.1765;0.9765 0.8471 0.1725;0.9725 0.851 0.1686;0.9686 0.8588 0.1647;0.9686 0.8667 0.1647;0.9647 0.8706 0.1608;0.9647 0.8784 0.1569;0.9608 0.8824 0.1569;0.9608 0.8902 0.1529;0.9608 0.898 0.149;0.9608 0.902 0.149;0.9608 0.9098 0.1451;0.9608 0.9137 0.1412;0.9608 0.9216 0.1373;0.9608 0.9255 0.1333;0.9608 0.9333 0.1294;0.9647 0.9373 0.1255;0.9647 0.9451 0.1216;0.9647 0.949 0.1176;0.9686 0.9569 0.1098;0.9686 0.9608 0.1059;0.9725 0.9686 0.102;0.9725 0.9725 0.0941;0.9765 0.9765 0.0863;0.9765 0.9843 0.0824];
            app.GelBoxUIFigure.Position = [100 100 1722 585];
            app.GelBoxUIFigure.Name = 'GelBox';
            app.GelBoxUIFigure.CloseRequestFcn = createCallbackFcn(app, @GelBoxUIFigureCloseRequest, true);

            % Create FileMenu
            app.FileMenu = uimenu(app.GelBoxUIFigure);
            app.FileMenu.Text = 'File';

            % Create LoadImageMenu
            app.LoadImageMenu = uimenu(app.FileMenu);
            app.LoadImageMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadImageButtonPushed, true);
            app.LoadImageMenu.Text = 'Load Image';

            % Create InvertImageMenu
            app.InvertImageMenu = uimenu(app.FileMenu);
            app.InvertImageMenu.MenuSelectedFcn = createCallbackFcn(app, @InvertImageButtonPushed, true);
            app.InvertImageMenu.Separator = 'on';
            app.InvertImageMenu.Text = 'Invert Image';

            % Create LoadAnalysisMenu
            app.LoadAnalysisMenu = uimenu(app.FileMenu);
            app.LoadAnalysisMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadAnalysisButtonPushed, true);
            app.LoadAnalysisMenu.Separator = 'on';
            app.LoadAnalysisMenu.Text = 'Load Analysis';

            % Create SaveAnalysisMenu
            app.SaveAnalysisMenu = uimenu(app.FileMenu);
            app.SaveAnalysisMenu.MenuSelectedFcn = createCallbackFcn(app, @SaveAnalysisButtonPushed, true);
            app.SaveAnalysisMenu.Separator = 'on';
            app.SaveAnalysisMenu.Text = 'Save Analysis';

            % Create ExportResultsMenu
            app.ExportResultsMenu = uimenu(app.FileMenu);
            app.ExportResultsMenu.MenuSelectedFcn = createCallbackFcn(app, @OutputButtonPushed, true);
            app.ExportResultsMenu.Separator = 'on';
            app.ExportResultsMenu.Text = 'Export Results';

            % Create DataAnalysisMenu
            app.DataAnalysisMenu = uimenu(app.GelBoxUIFigure);
            app.DataAnalysisMenu.Enable = 'off';
            app.DataAnalysisMenu.Text = 'Data Analysis';

            % Create GelImageFileInformationMenu
            app.GelImageFileInformationMenu = uimenu(app.DataAnalysisMenu);
            app.GelImageFileInformationMenu.MenuSelectedFcn = createCallbackFcn(app, @GelImageFileInformationMenuSelected, true);
            app.GelImageFileInformationMenu.Text = 'Gel Image File Information';

            % Create SelectedBoxInformationMenu
            app.SelectedBoxInformationMenu = uimenu(app.DataAnalysisMenu);
            app.SelectedBoxInformationMenu.MenuSelectedFcn = createCallbackFcn(app, @SelectedBoxInformationMenuSelected, true);
            app.SelectedBoxInformationMenu.Separator = 'on';
            app.SelectedBoxInformationMenu.Text = 'Selected Box Information';

            % Create SummaryPlotMenu
            app.SummaryPlotMenu = uimenu(app.DataAnalysisMenu);
            app.SummaryPlotMenu.MenuSelectedFcn = createCallbackFcn(app, @AnalysisSummaryPlotMenuSelected, true);
            app.SummaryPlotMenu.Separator = 'on';
            app.SummaryPlotMenu.Text = 'Summary Plot';

            % Create OpticalDensitiesPanel
            app.OpticalDensitiesPanel = uipanel(app.GelBoxUIFigure);
            app.OpticalDensitiesPanel.Title = 'Optical Densities';
            app.OpticalDensitiesPanel.Position = [875 299 836 278];

            % Create raw_density
            app.raw_density = uiaxes(app.OpticalDensitiesPanel);
            xlabel(app.raw_density, 'Optical Density (A.U.)')
            ylabel(app.raw_density, 'Pixel')
            app.raw_density.Position = [163 11 254 209];

            % Create background_corrected_raw_density
            app.background_corrected_raw_density = uiaxes(app.OpticalDensitiesPanel);
            xlabel(app.background_corrected_raw_density, 'Optical Density (A.U.)')
            ylabel(app.background_corrected_raw_density, 'Pixel')
            app.background_corrected_raw_density.YColor = [0.9412 0.9412 0.9412];
            app.background_corrected_raw_density.Position = [401 11 254 209];

            % Create box_inset
            app.box_inset = uiaxes(app.OpticalDensitiesPanel);
            xlabel(app.box_inset, 'Optical Density')
            ylabel(app.box_inset, 'Pixel')
            app.box_inset.XColor = [0.9412 0.9412 0.9412];
            app.box_inset.YColor = [0.9412 0.9412 0.9412];
            app.box_inset.Position = [4 11 150 209];

            % Create BoxZoomLabel
            app.BoxZoomLabel = uilabel(app.OpticalDensitiesPanel);
            app.BoxZoomLabel.HorizontalAlignment = 'center';
            app.BoxZoomLabel.WordWrap = 'on';
            app.BoxZoomLabel.Position = [34 224 126 22];
            app.BoxZoomLabel.Text = 'Box Zoom';

            % Create RawOpticalDensityLabel_2
            app.RawOpticalDensityLabel_2 = uilabel(app.OpticalDensitiesPanel);
            app.RawOpticalDensityLabel_2.HorizontalAlignment = 'center';
            app.RawOpticalDensityLabel_2.WordWrap = 'on';
            app.RawOpticalDensityLabel_2.Position = [238 224 126 22];
            app.RawOpticalDensityLabel_2.Text = 'Raw Optical Density';

            % Create BackgroundCorrectedOpticalDensityLabel_2
            app.BackgroundCorrectedOpticalDensityLabel_2 = uilabel(app.OpticalDensitiesPanel);
            app.BackgroundCorrectedOpticalDensityLabel_2.HorizontalAlignment = 'center';
            app.BackgroundCorrectedOpticalDensityLabel_2.WordWrap = 'on';
            app.BackgroundCorrectedOpticalDensityLabel_2.Position = [477 221 126 28];
            app.BackgroundCorrectedOpticalDensityLabel_2.Text = 'Background Corrected Optical Density';

            % Create TotalAreaEditFieldLabel
            app.TotalAreaEditFieldLabel = uilabel(app.OpticalDensitiesPanel);
            app.TotalAreaEditFieldLabel.Position = [661 146 59 22];
            app.TotalAreaEditFieldLabel.Text = 'Total Area';

            % Create TotalAreaField
            app.TotalAreaField = uieditfield(app.OpticalDensitiesPanel, 'numeric');
            app.TotalAreaField.ValueDisplayFormat = '%.2f';
            app.TotalAreaField.Editable = 'off';
            app.TotalAreaField.HorizontalAlignment = 'center';
            app.TotalAreaField.Position = [732 146 85 22];

            % Create BackgroundAreaLabel
            app.BackgroundAreaLabel = uilabel(app.OpticalDensitiesPanel);
            app.BackgroundAreaLabel.WordWrap = 'on';
            app.BackgroundAreaLabel.Position = [661 106 81 30];
            app.BackgroundAreaLabel.Text = 'Background Area';

            % Create BackgroundAreaField
            app.BackgroundAreaField = uieditfield(app.OpticalDensitiesPanel, 'numeric');
            app.BackgroundAreaField.ValueDisplayFormat = '%.2f';
            app.BackgroundAreaField.Editable = 'off';
            app.BackgroundAreaField.HorizontalAlignment = 'center';
            app.BackgroundAreaField.Position = [732 110 85 22];

            % Create GelImagePanel
            app.GelImagePanel = uipanel(app.GelBoxUIFigure);
            app.GelImagePanel.Title = 'Gel Image';
            app.GelImagePanel.Position = [6 10 860 567];

            % Create gel_image_axes
            app.gel_image_axes = uiaxes(app.GelImagePanel);
            app.gel_image_axes.XTick = [];
            app.gel_image_axes.YTick = [];
            app.gel_image_axes.Box = 'on';
            app.gel_image_axes.Position = [12 13 837 487];

            % Create DeleteBoxButton
            app.DeleteBoxButton = uibutton(app.GelImagePanel, 'push');
            app.DeleteBoxButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteBoxButtonPushed, true);
            app.DeleteBoxButton.Enable = 'off';
            app.DeleteBoxButton.Position = [127 514 101 22];
            app.DeleteBoxButton.Text = 'Delete Box';

            % Create NewBoxButton
            app.NewBoxButton = uibutton(app.GelImagePanel, 'push');
            app.NewBoxButton.ButtonPushedFcn = createCallbackFcn(app, @NewBoxButtonPushed, true);
            app.NewBoxButton.Position = [13 514 100 22];
            app.NewBoxButton.Text = 'New Box';

            % Create BoxSelectionDropDownLabel
            app.BoxSelectionDropDownLabel = uilabel(app.GelImagePanel);
            app.BoxSelectionDropDownLabel.HorizontalAlignment = 'center';
            app.BoxSelectionDropDownLabel.Position = [231 515 98 22];
            app.BoxSelectionDropDownLabel.Text = 'Box Selection';

            % Create BoxSelectionDropDown
            app.BoxSelectionDropDown = uidropdown(app.GelImagePanel);
            app.BoxSelectionDropDown.Items = {};
            app.BoxSelectionDropDown.ValueChangedFcn = createCallbackFcn(app, @BoxSelectionDropDownValueChanged, true);
            app.BoxSelectionDropDown.Placeholder = 'No data';
            app.BoxSelectionDropDown.Position = [329 515 100 22];
            app.BoxSelectionDropDown.Value = {};

            % Create FittingPanel
            app.FittingPanel = uipanel(app.GelBoxUIFigure);
            app.FittingPanel.Title = 'Fitting';
            app.FittingPanel.Position = [875 10 836 277];

            % Create background_corrected_raw_density_fit
            app.background_corrected_raw_density_fit = uiaxes(app.FittingPanel);
            xlabel(app.background_corrected_raw_density_fit, 'Optical Density (A.U.)')
            ylabel(app.background_corrected_raw_density_fit, 'Pixel')
            app.background_corrected_raw_density_fit.YColor = [0.9412 0.9412 0.9412];
            app.background_corrected_raw_density_fit.Position = [251 7 254 209];

            % Create raw_density_fit
            app.raw_density_fit = uiaxes(app.FittingPanel);
            xlabel(app.raw_density_fit, 'Optical Density (A.U.)')
            ylabel(app.raw_density_fit, 'Pixel')
            app.raw_density_fit.Position = [17 7 254 209];

            % Create RawOpticalDensityLabel
            app.RawOpticalDensityLabel = uilabel(app.FittingPanel);
            app.RawOpticalDensityLabel.HorizontalAlignment = 'center';
            app.RawOpticalDensityLabel.WordWrap = 'on';
            app.RawOpticalDensityLabel.Position = [82 224 126 22];
            app.RawOpticalDensityLabel.Text = 'Raw Optical Density';

            % Create BackgroundCorrectedOpticalDensityLabel
            app.BackgroundCorrectedOpticalDensityLabel = uilabel(app.FittingPanel);
            app.BackgroundCorrectedOpticalDensityLabel.HorizontalAlignment = 'center';
            app.BackgroundCorrectedOpticalDensityLabel.WordWrap = 'on';
            app.BackgroundCorrectedOpticalDensityLabel.Position = [331 221 126 28];
            app.BackgroundCorrectedOpticalDensityLabel.Text = 'Background Corrected Optical Density';

            % Create DrawFittingCheckBox
            app.DrawFittingCheckBox = uicheckbox(app.FittingPanel);
            app.DrawFittingCheckBox.Text = 'Draw Fitting';
            app.DrawFittingCheckBox.Position = [672 154 86 22];

            % Create RsquaredLabel
            app.RsquaredLabel = uilabel(app.FittingPanel);
            app.RsquaredLabel.WordWrap = 'on';
            app.RsquaredLabel.Position = [512 154 72 22];
            app.RsquaredLabel.Text = 'R - squared';

            % Create rsquaredField
            app.rsquaredField = uieditfield(app.FittingPanel, 'numeric');
            app.rsquaredField.ValueDisplayFormat = '%.3f';
            app.rsquaredField.Editable = 'off';
            app.rsquaredField.HorizontalAlignment = 'center';
            app.rsquaredField.Position = [588 154 70 22];

            % Create BandAreaLabel_1
            app.BandAreaLabel_1 = uilabel(app.FittingPanel);
            app.BandAreaLabel_1.HorizontalAlignment = 'center';
            app.BandAreaLabel_1.WordWrap = 'on';
            app.BandAreaLabel_1.Position = [512 112 72 28];
            app.BandAreaLabel_1.Text = 'Band Area 1 (Red)';

            % Create BandArea_1
            app.BandArea_1 = uieditfield(app.FittingPanel, 'numeric');
            app.BandArea_1.ValueDisplayFormat = '%.2f';
            app.BandArea_1.Editable = 'off';
            app.BandArea_1.HorizontalAlignment = 'center';
            app.BandArea_1.Position = [587 116 70 22];

            % Create BandRelativeAreaLabel_1
            app.BandRelativeAreaLabel_1 = uilabel(app.FittingPanel);
            app.BandRelativeAreaLabel_1.WordWrap = 'on';
            app.BandRelativeAreaLabel_1.Enable = 'off';
            app.BandRelativeAreaLabel_1.Position = [671 105 81 42];
            app.BandRelativeAreaLabel_1.Text = 'Band Relative Area 1 (Red)';

            % Create BandRelativeArea_1
            app.BandRelativeArea_1 = uieditfield(app.FittingPanel, 'numeric');
            app.BandRelativeArea_1.ValueDisplayFormat = '%.2f';
            app.BandRelativeArea_1.Editable = 'off';
            app.BandRelativeArea_1.HorizontalAlignment = 'center';
            app.BandRelativeArea_1.Enable = 'off';
            app.BandRelativeArea_1.Position = [756 115 70 22];

            % Create BandAreaLabel_2
            app.BandAreaLabel_2 = uilabel(app.FittingPanel);
            app.BandAreaLabel_2.HorizontalAlignment = 'center';
            app.BandAreaLabel_2.WordWrap = 'on';
            app.BandAreaLabel_2.Enable = 'off';
            app.BandAreaLabel_2.Position = [511 73 72 28];
            app.BandAreaLabel_2.Text = 'Band Area 2 (Blue)';

            % Create BandArea_2
            app.BandArea_2 = uieditfield(app.FittingPanel, 'numeric');
            app.BandArea_2.ValueDisplayFormat = '%.2f';
            app.BandArea_2.Editable = 'off';
            app.BandArea_2.HorizontalAlignment = 'center';
            app.BandArea_2.Enable = 'off';
            app.BandArea_2.Position = [587 76 70 22];

            % Create BandRelativeAreaLabel_2
            app.BandRelativeAreaLabel_2 = uilabel(app.FittingPanel);
            app.BandRelativeAreaLabel_2.WordWrap = 'on';
            app.BandRelativeAreaLabel_2.Enable = 'off';
            app.BandRelativeAreaLabel_2.Position = [670 66 81 42];
            app.BandRelativeAreaLabel_2.Text = 'Band Relative Area 2 (Blue)';

            % Create BandRelativeArea_2
            app.BandRelativeArea_2 = uieditfield(app.FittingPanel, 'numeric');
            app.BandRelativeArea_2.ValueDisplayFormat = '%.2f';
            app.BandRelativeArea_2.Editable = 'off';
            app.BandRelativeArea_2.HorizontalAlignment = 'center';
            app.BandRelativeArea_2.Enable = 'off';
            app.BandRelativeArea_2.Position = [756 76 70 22];

            % Create BandAreaLabel_3
            app.BandAreaLabel_3 = uilabel(app.FittingPanel);
            app.BandAreaLabel_3.HorizontalAlignment = 'center';
            app.BandAreaLabel_3.WordWrap = 'on';
            app.BandAreaLabel_3.Enable = 'off';
            app.BandAreaLabel_3.Position = [511 35 72 30];
            app.BandAreaLabel_3.Text = 'Band Area 3 (Yellow)';

            % Create BandArea_3
            app.BandArea_3 = uieditfield(app.FittingPanel, 'numeric');
            app.BandArea_3.ValueDisplayFormat = '%.2f';
            app.BandArea_3.Editable = 'off';
            app.BandArea_3.HorizontalAlignment = 'center';
            app.BandArea_3.Enable = 'off';
            app.BandArea_3.Position = [587 40 70 22];

            % Create BandRelativeAreaLabel_3
            app.BandRelativeAreaLabel_3 = uilabel(app.FittingPanel);
            app.BandRelativeAreaLabel_3.WordWrap = 'on';
            app.BandRelativeAreaLabel_3.Enable = 'off';
            app.BandRelativeAreaLabel_3.Position = [668 29 81 42];
            app.BandRelativeAreaLabel_3.Text = 'Band Relative Area 3 (Yellow)';

            % Create BandRelativeArea_3
            app.BandRelativeArea_3 = uieditfield(app.FittingPanel, 'numeric');
            app.BandRelativeArea_3.ValueDisplayFormat = '%.2f';
            app.BandRelativeArea_3.Editable = 'off';
            app.BandRelativeArea_3.HorizontalAlignment = 'center';
            app.BandRelativeArea_3.Enable = 'off';
            app.BandRelativeArea_3.Position = [756 39 70 22];

            % Create NumberofBandsDropDownLabel
            app.NumberofBandsDropDownLabel = uilabel(app.FittingPanel);
            app.NumberofBandsDropDownLabel.Position = [512 194 99 22];
            app.NumberofBandsDropDownLabel.Text = 'Number of Bands';

            % Create NumberofBandsDropDown
            app.NumberofBandsDropDown = uidropdown(app.FittingPanel);
            app.NumberofBandsDropDown.Items = {'1', '2', '3'};
            app.NumberofBandsDropDown.ValueChangedFcn = createCallbackFcn(app, @NumberofBandsDropDownValueChanged, true);
            app.NumberofBandsDropDown.Position = [614 194 42 22];
            app.NumberofBandsDropDown.Value = '1';

            % Create FittingParametersButton
            app.FittingParametersButton = uibutton(app.FittingPanel, 'push');
            app.FittingParametersButton.ButtonPushedFcn = createCallbackFcn(app, @FittingParametersButtonPushed, true);
            app.FittingParametersButton.Position = [666 193 113 23];
            app.FittingParametersButton.Text = 'Fitting Parameters';

            % Show the figure after all components are created
            app.GelBoxUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GelBox_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.GelBoxUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.GelBoxUIFigure)
        end
    end
end