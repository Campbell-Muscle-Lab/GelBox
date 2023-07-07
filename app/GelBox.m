
classdef GelBox < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GelBoxUIFigure               matlab.ui.Figure
        FileMenu                     matlab.ui.container.Menu
        LoadImageMenu                matlab.ui.container.Menu
        InvertImageMenu              matlab.ui.container.Menu
        LoadAnalysisMenu             matlab.ui.container.Menu
        SaveAnalysisMenu             matlab.ui.container.Menu
        OutputMenu                   matlab.ui.container.Menu
        GelImageFileInformationMenu  matlab.ui.container.Menu
        SelectedBoxInformationMenu   matlab.ui.container.Menu
        AnalysisSummaryPlotMenu      matlab.ui.container.Menu
        FittingPanel                 matlab.ui.container.Panel
        FittingOptionsButton         matlab.ui.control.Button
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

                    switch num_of_bands
                        case 1
                            app.BandRelativeArea_1.Enable = 0;
                            [x_bands,x_fit,r_squared] = ...
                            FitGaussian(app,y,x,x_back,num_of_bands);
                        case 2
                            if ~app.BandRelativeArea_1.Enable
                                app.BandRelativeArea_1.Enable = 1;
                            end
                            [x_bands,x_fit,r_squared] = ...
                            Fit2Gaussian(app,y,x,x_back,num_of_bands);
                        case 3
                            if ~app.BandRelativeArea_1.Enable
                                app.BandRelativeArea_1.Enable = 1;
                            end
                            [x_bands,x_fit,r_squared] = ...
                            Fit3Gaussian(app,y,x,x_back,num_of_bands);
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
                       legend(app.raw_density,'','Baseline', ...
                           'Location','best')


                       xticks(app.raw_density_fit,x_ticks);
                       app.raw_density_fit.XAxis.Exponent = 0;

                       xlim(app.raw_density_fit,[0 x_t_end]);
                       ylim(app.raw_density_fit,[1 max(y)]);
                       legend(app.raw_density)


                       cla(app.background_corrected_raw_density)
                       plot(app.background_corrected_raw_density, ...
                           x-x_back',y,'-.k',"LineWidth",2)
                       hold(app.background_corrected_raw_density,"on")
                       plot(app.background_corrected_raw_density, ...
                           zeros(1,numel(y)),y,'-.m',"LineWidth",2)
                       legend(app.background_corrected_raw_density,'','Baseline', ...
                           'Location','best')
                       ylim(app.background_corrected_raw_density, ...
                           [1 max(y)]);
                       
                       x_pow = ceil(log10(max(x-x_back')));
                       x_tick_rounder = 10^(x_pow - 1);
                       x_t_end = ceil(max(x-x_back')/x_tick_rounder)*x_tick_rounder;
                       if min(x-x_back') < 0
                           x_tick_rounder = -10^(x_pow - 1)*0.75;
                       else
                           x_tick_rounder = 10^(x_pow - 1)*0.75;
                       end
                       x_t_beginning = ceil(min(x-x_back')/x_tick_rounder)*x_tick_rounder;
                       x_t_mid = round(x_t_end/2);
                       if x_t_beginning < 0
                           x_ticks = [x_t_beginning 0 x_t_mid x_t_end];
                       else
                           x_ticks = [x_t_beginning x_t_mid x_t_end];
                       end
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
                               color = {'r','b','g'};
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
                       plot(app.raw_density_fit, ...
                           x,y,"Color",'k',"LineWidth",2)
                       plot(app.background_corrected_raw_density_fit, ...
                           x-x_back',y,'-.k',"LineWidth",2)
                       for j = 1 : num_of_bands
                            patch(app.raw_density_fit, ...
                                x_back+x_bands(j,:), ...
                                y,color{j},'FaceAlpha',0.25, ...
                                'EdgeColor',color{j},'EdgeAlpha',0.25, ...
                                'LineWidth',2)
                            patch(app.background_corrected_raw_density_fit, ...
                                x_bands(j,:), ...
                                y,color{j},'FaceAlpha',0.25, ...
                                'EdgeColor',color{j},'EdgeAlpha',0.25, ...
                                'LineWidth',2)
                       end
                       ylim(app.raw_density_fit, ...
                           [1 max(y)]);
                       ylim(app.background_corrected_raw_density_fit, ...
                           [1 max(y)]);
                    end
                end
                app.gel_data.d_box = d;
            end
        end
        
       function [y_bands, y_fit,r_squared] = FitGaussian(app,x,y,y_back, ...
                no_of_bands)

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

            par = [first_curve_x_estimate ...
                   first_curve_shape_estimate ...
                   first_curve_amp_estimate ...
                   first_curve_skew_estimate ...
                   ];

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

                if any(par<0)
                    trial_e = 10^6;
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
                y_back,no_of_bands)
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

            
            par = [first_curve_x_estimate ...
                   first_curve_shape_estimate ...
                   first_curve_amp_estimate ...
                   first_curve_skew_estimate ...
                   second_curve_x_estimate ...
                   second_curve_amp_estimate ...
                   ];

            if ~app.fitting_options.shared_shape && ~app.fitting_options.shared_skewness
                par(7) = second_curve_shape_estimate;
                par(8) = second_curve_skew_estimate;
            elseif ~app.fitting_options.shared_shape && app.fitting_options.shared_skewness
                par(7) = second_curve_shape_estimate;
            elseif app.fitting_options.shared_shape && ~app.fitting_options.shared_skewness
                par(7) = second_curve_skew_estimate;
            end


            j = 1;
            e = [];

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
                        trial_e = 10^12;
                    end
                end

                if any(par<0)
                    trial_e = 10^12;
                end
                
                for i = 1:2
                    areas(i) = trapz(x,y_bands(i,:));
                end

                if any(areas<0)
                    trial_e = 10^12;
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
                                
                if ~app.fitting_options.shared_shape && ...
                    ~app.fitting_options.shared_skewness
                    curve_shape2 = par(7);
                    skew2 = par(8);
                elseif ~app.fitting_options.shared_shape && ...
                    app.fitting_options.shared_skewness
                    curve_shape2 = par(7);
                    skew2 = skew1;
                elseif app.fitting_options.shared_shape && ...
                    ~app.fitting_options.shared_skewness
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

        function [y_bands, y_fit,r_squared] = Fit3Gaussian(app,x,y,y_back, ...
                no_of_bands)
            peaks=find_peaks('x',x, ...
                             'y',y, ...
                             'min_rel_delta_y',0.05, ...
                             'min_x_index_spacing',2);
            
            if numel(peaks.max_indices) == no_of_bands
                first_curve_x_estimate=peaks.max_indices(1);
                second_curve_x_estimate=peaks.max_indices(2);
                third_curve_x_estimate=peaks.max_indices(3);
            else
                first_curve_x_estimate = 0.3*length(x);
                second_curve_x_estimate = 0.6*length(x);
                third_curve_x_estimate = 0.9*length(x);
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

            par = [first_curve_x_estimate ...
                   first_curve_shape_estimate ...
                   first_curve_amp_estimate ...
                   first_curve_skew_estimate ...
                   second_curve_x_estimate ...
                   second_curve_amp_estimate ...
                   third_curve_x_estimate ...
                   third_curve_amp_estimate ...
                   ];

            if ~app.fitting_options.shared_shape && ~app.fitting_options.shared_skewness
                par(9) = second_curve_shape_estimate;
                par(10) = second_curve_skew_estimate;
                par(11) = third_curve_shape_estimate;
                par(12) = third_curve_skew_estimate;
            elseif ~app.fitting_options.shared_shape && app.fitting_options.shared_skewness
                par(9) = second_curve_shape_estimate;
                par(10) = third_curve_shape_estimate;
            elseif app.fitting_options.shared_shape && ~app.fitting_options.shared_skewness
                par(9) = second_curve_skew_estimate;
                par(10) = third_curve_skew_estimate;

            end


            j = 1;
            e = [];

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
                        trial_e = 10^12;
                    end
                end

                if any(par<0)
                    trial_e = 10^12;
                end
                
                for i = 1:3
                    areas(i) = trapz(x,y_bands(i,:));
                end

                if any(areas<0)
                    trial_e = 10^12;
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
                
                if ~app.fitting_options.shared_shape && ...
                    ~app.fitting_options.shared_skewness
                    curve_shape2 = par(9);
                    skew2 = par(10);
                    curve_shape3 = par(11);
                    skew3 = par(12);
                elseif ~app.fitting_options.shared_shape && ...
                    app.fitting_options.shared_skewness
                    curve_shape2 = par(9);
                    skew2 = skew1;
                    curve_shape3 = par(10);
                    skew3 = skew1;
                elseif app.fitting_options.shared_shape && ...
                    ~app.fitting_options.shared_skewness
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
        
        function UpdateFittingOptions(app,shape,skewness)
            app.fitting_options.shared_shape = shape;
            app.fitting_options.shared_skewness = skewness;
            UpdateDisplay(app)
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
            app.NumberofBandsDropDownValueChanged;

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
            writelines(evalc('type(mfilename(''fullpath'')+".mlapp")'),mfilename('fullpath')+".m");
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
                app.GelImageFileInformationMenu.Enable = 1;
                app.SelectedBoxInformationMenu.Enable = 1;
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
                
                app.DeleteBoxButton.Enable = 1;
                app.OutputButton.Enable = 1;
                app.GelImageFileInformationMenu.Enable = 1;
                app.SelectedBoxInformationMenu.Enable = 1;

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
                control_strings = [];

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

            % Nested function
            function new_box_position2(evt);
                app.gel_data = guidata(gui.Window);
                if (isfield(app.gel_data,'box_position'))
                    box_position = app.gel_data.box_position;
                    [r,c]=size(box_position);
                    if (r>=n)&(~isequal(box_position(n,:),evt.CurrentPosition))
                        update_display(gui,n);
                    end
                else
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
            
            app.DeleteBoxButton.Enable = 1;
            app.SelectedBoxInformationMenu.Enable = 1;
            UpdateDisplay(app);

            function new_box_position(evt);
                if (isfield(app.gel_data,'box_position'))
                    box_position = app.gel_data.box_position;
                    [r,c]=size(box_position);
                    if (r>=n)&(~isequal(box_position(n,:),evt.CurrentPosition))
                        UpdateDisplay(app);
                    end
                else
                        UpdateDisplay(app);
                end
            end

        end

        % Value changed function: NumberofBandsDropDown
        function NumberofBandsDropDownValueChanged(app, event)
            value = app.NumberofBandsDropDown.Value;
            
            switch value
                case '1'
                    try
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

        % Menu selected function: OutputMenu
        function OutputButtonPushed(app, event)
            d = [];
            d.image_file{1} = app.gel_data.image_file_string;
            n = numel(app.gel_data.box_handle);

            for i=1:n
                d.band(i) = i;
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
                d.num_of_bands(i) = band_orientation;
                d.r_squared(i) = app.gel_data.summary(i).r_squared;
            end

            [file_string,path_string] = uiputfile2( ...
                {'*.xlsx','Excel file'},'Select file for results');
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

        % Button pushed function: FittingOptionsButton
        function FittingOptionsButtonPushed(app, event)
            app.FittingOptions = FittingOptionsDialog(app);

        end

        % Menu selected function: AnalysisSummaryPlotMenu
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
            colormap(app.GelBoxUIFigure, 'parula');
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
            app.InvertImageMenu.Text = 'Invert Image';

            % Create LoadAnalysisMenu
            app.LoadAnalysisMenu = uimenu(app.FileMenu);
            app.LoadAnalysisMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadAnalysisButtonPushed, true);
            app.LoadAnalysisMenu.Text = 'Load Analysis';

            % Create SaveAnalysisMenu
            app.SaveAnalysisMenu = uimenu(app.FileMenu);
            app.SaveAnalysisMenu.MenuSelectedFcn = createCallbackFcn(app, @SaveAnalysisButtonPushed, true);
            app.SaveAnalysisMenu.Text = 'Save Analysis';

            % Create OutputMenu
            app.OutputMenu = uimenu(app.FileMenu);
            app.OutputMenu.MenuSelectedFcn = createCallbackFcn(app, @OutputButtonPushed, true);
            app.OutputMenu.Text = 'Output';

            % Create GelImageFileInformationMenu
            app.GelImageFileInformationMenu = uimenu(app.GelBoxUIFigure);
            app.GelImageFileInformationMenu.MenuSelectedFcn = createCallbackFcn(app, @GelImageFileInformationMenuSelected, true);
            app.GelImageFileInformationMenu.Enable = 'off';
            app.GelImageFileInformationMenu.Text = 'Gel Image File Information';

            % Create SelectedBoxInformationMenu
            app.SelectedBoxInformationMenu = uimenu(app.GelBoxUIFigure);
            app.SelectedBoxInformationMenu.MenuSelectedFcn = createCallbackFcn(app, @SelectedBoxInformationMenuSelected, true);
            app.SelectedBoxInformationMenu.Enable = 'off';
            app.SelectedBoxInformationMenu.Text = 'Selected Box Information';

            % Create AnalysisSummaryPlotMenu
            app.AnalysisSummaryPlotMenu = uimenu(app.GelBoxUIFigure);
            app.AnalysisSummaryPlotMenu.MenuSelectedFcn = createCallbackFcn(app, @AnalysisSummaryPlotMenuSelected, true);
            app.AnalysisSummaryPlotMenu.Text = 'Analysis Summary Plot';

            % Create OpticalDensitiesPanel
            app.OpticalDensitiesPanel = uipanel(app.GelBoxUIFigure);
            app.OpticalDensitiesPanel.Title = 'Optical Densities';
            app.OpticalDensitiesPanel.Position = [7 299 836 278];

            % Create raw_density
            app.raw_density = uiaxes(app.OpticalDensitiesPanel);
            xlabel(app.raw_density, 'Optical Density')
            ylabel(app.raw_density, 'Pixel')
            app.raw_density.Position = [163 11 254 209];

            % Create background_corrected_raw_density
            app.background_corrected_raw_density = uiaxes(app.OpticalDensitiesPanel);
            xlabel(app.background_corrected_raw_density, 'Optical Density')
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
            app.GelImagePanel.Position = [853 10 860 567];

            % Create gel_image_axes
            app.gel_image_axes = uiaxes(app.GelImagePanel);
            app.gel_image_axes.XTick = [];
            app.gel_image_axes.YTick = [];
            app.gel_image_axes.Box = 'on';
            app.gel_image_axes.Visible = 'off';
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
            app.FittingPanel.Position = [7 10 836 277];

            % Create background_corrected_raw_density_fit
            app.background_corrected_raw_density_fit = uiaxes(app.FittingPanel);
            xlabel(app.background_corrected_raw_density_fit, 'Optical Density')
            ylabel(app.background_corrected_raw_density_fit, 'Pixel')
            app.background_corrected_raw_density_fit.YColor = [0.9412 0.9412 0.9412];
            app.background_corrected_raw_density_fit.Position = [251 7 254 209];

            % Create raw_density_fit
            app.raw_density_fit = uiaxes(app.FittingPanel);
            xlabel(app.raw_density_fit, 'Optical Density')
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
            app.rsquaredField.ValueDisplayFormat = '%.2f';
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
            app.BandAreaLabel_3.Position = [511 37 72 28];
            app.BandAreaLabel_3.Text = 'Band Area 3 (Green)';

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
            app.BandRelativeAreaLabel_3.Text = 'Band Relative Area 3 (Green)';

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

            % Create FittingOptionsButton
            app.FittingOptionsButton = uibutton(app.FittingPanel, 'push');
            app.FittingOptionsButton.ButtonPushedFcn = createCallbackFcn(app, @FittingOptionsButtonPushed, true);
            app.FittingOptionsButton.Position = [672 194 100 22];
            app.FittingOptionsButton.Text = 'Fitting Options';

            % Show the figure after all components are created
            app.GelBoxUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GelBox

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

