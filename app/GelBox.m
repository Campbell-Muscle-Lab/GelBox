
classdef GelBox < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GelBoxUIFigure                matlab.ui.Figure
        SelectedBoxFittingPanel       matlab.ui.container.Panel
        BandRelativeArea_3            matlab.ui.control.NumericEditField
        BandRelativeAreaLabel_3       matlab.ui.control.Label
        BandArea_3                    matlab.ui.control.NumericEditField
        BandAreaLabel_3               matlab.ui.control.Label
        BandRelativeArea_2            matlab.ui.control.NumericEditField
        BandRelativeAreaLabel_2       matlab.ui.control.Label
        BandArea_2                    matlab.ui.control.NumericEditField
        BandAreaLabel_2               matlab.ui.control.Label
        BandRelativeArea_1            matlab.ui.control.NumericEditField
        BandRelativeAreaLabel_1       matlab.ui.control.Label
        BandArea_1                    matlab.ui.control.NumericEditField
        BandAreaLabel_1               matlab.ui.control.Label
        rsquaredField                 matlab.ui.control.NumericEditField
        RsquaredLabel                 matlab.ui.control.Label
        DrawFittingCheckBox           matlab.ui.control.CheckBox
        BackgroundCorrectedOpticalDensityLabel  matlab.ui.control.Label
        RawOpticalDensityLabel        matlab.ui.control.Label
        raw_density_fit               matlab.ui.control.UIAxes
        background_corrected_raw_density_fit  matlab.ui.control.UIAxes
        GelImagePanel                 matlab.ui.container.Panel
        gel_image_axes                matlab.ui.control.UIAxes
        GelImageFileInformationPanel  matlab.ui.container.Panel
        ImFInfoArea                   matlab.ui.control.TextArea
        SelectedBoxOpticalDensitiesPanel  matlab.ui.container.Panel
        BackgroundAreaField           matlab.ui.control.NumericEditField
        BackgroundAreaLabel           matlab.ui.control.Label
        TotalAreaField                matlab.ui.control.NumericEditField
        TotalAreaEditFieldLabel       matlab.ui.control.Label
        BackgroundCorrectedOpticalDensityLabel_2  matlab.ui.control.Label
        RawOpticalDensityLabel_2      matlab.ui.control.Label
        BoxZoomLabel                  matlab.ui.control.Label
        box_inset                     matlab.ui.control.UIAxes
        background_corrected_raw_density  matlab.ui.control.UIAxes
        raw_density                   matlab.ui.control.UIAxes
        SelectedBoxInformationPanel   matlab.ui.container.Panel
        SelectedBoxInformationArea    matlab.ui.control.TextArea
        AnalysisControlsPanel         matlab.ui.container.Panel
        BoxSelectionDropDown          matlab.ui.control.DropDown
        BoxSelectionDropDownLabel     matlab.ui.control.Label
        DeleteBoxButton               matlab.ui.control.Button
        NewBoxButton                  matlab.ui.control.Button
        NumberofBandsDropDown         matlab.ui.control.DropDown
        NumberofBandsDropDownLabel    matlab.ui.control.Label
        FileControlsPanel             matlab.ui.container.Panel
        OutputButton                  matlab.ui.control.Button
        SaveAnalysisButton            matlab.ui.control.Button
        LoadAnalysisButton            matlab.ui.control.Button
        InvertImageButton             matlab.ui.control.Button
        LoadImageButton               matlab.ui.control.Button
    end

    
    properties (Access = public)
        gel_data % Description
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
                       legend(app.raw_density,'',{'Baseline'}, ...
                           'Location','northwest')


                       xticks(app.raw_density_fit,x_ticks);
                       app.raw_density_fit.XAxis.Exponent = 0;

                       xlim(app.raw_density_fit,[0 x_t_end]);
                       ylim(app.raw_density_fit,[1 max(y)]);
                       legend(app.raw_density)



                       plot(app.background_corrected_raw_density, ...
                           x-x_back',y,'-.k',"LineWidth",2)
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
                            
                       app.SelectedBoxInformationArea.Value = printstruct(d.box(i));

                    end
                end
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

            par = [first_curve_x_estimate ...
                   first_curve_shape_estimate ...
                   first_curve_amp_estimate ...
                   first_curve_skew_estimate ...
                   second_curve_x_estimate ...
                   second_curve_amp_estimate ...
                   second_curve_shape_estimate ...
                   ];

            j = 1;
            e = [];

            opts=optimset('fminsearch');
            opts.Display='off';
            opts.MaxIter=1000;
            opts.MaxFunEvals=10000;

            [p_result,fval,exitflag,output] = fminsearch(@profile_error_2gaussian, par, opts);

            r_squared = calculate_r_squared(y',y_fit+y_back);

            function trial_e = profile_error_2gaussian(par)

                [y_bands,y_fit] = calculate_2profile(x,par);

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
            function [y_bands,y_fit] = calculate_2profile(x,par)
                
                x1 = par(1);
                curve_shape1 = par(2);
                amp1 = par(3);
                skew1 = par(4);
                x2 = par(5);
                amp2 = par(6);
                curve_shape2 = par(7);
                
                %         x1
                y_first = skewed_Gaussian(x,x1,curve_shape1,amp1,skew1);
                y_second = skewed_Gaussian(x,x2,curve_shape2,amp2,skew1);

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

            third_curve_amp_estimate = first_curve_amp_estimate;
            third_curve_shape_estimate = alfa_estimate;

            par = [first_curve_x_estimate ...
                   first_curve_shape_estimate ...
                   first_curve_amp_estimate ...
                   first_curve_skew_estimate ...
                   second_curve_x_estimate ...
                   second_curve_amp_estimate ...
                   third_curve_x_estimate ...
                   third_curve_amp_estimate ...
                   ];

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
            function [y_bands,y_fit] = calculate_3profile(x,par)
                
                x1 = par(1);
                curve_shape1 = par(2);
                amp1 = par(3);
                skew1 = par(4);
                x2 = par(5);
                amp2 = par(6);
                x3 = par(7);
                amp3 = par(8);
                
                %         x1
                y_first = skewed_Gaussian(x,x1,curve_shape1,amp1,skew1);
                y_second = skewed_Gaussian(x,x2,curve_shape1,amp2,skew1);
                y_third = skewed_Gaussian(x,x3,curve_shape1,amp3,skew1);


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
    end
    
    methods (Access = private)
        
        function ResetDisplay(app)
            % Reset the file controls
            app.SaveAnalysisButton.Enable = 0;
            app.OutputButton.Enable = 0;

            % Clear the figure displays
            cla(app.box_inset)
            cla(app.raw_density)
            cla(app.background_corrected_raw_density)
            cla(app.raw_density_fit)
            cla(app.background_corrected_raw_density_fit)

            % Clear the text displays
            app.SelectedBoxInformationArea.Value = '';
            app.ImFInfoArea.Value = '';

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

        end

        % Button pushed function: LoadImageButton
        function LoadImageButtonPushed(app, event)
            [file_string,path_string]=uigetfile2( ...
                {'*.png','PNG';'*.tif','TIF'}, ...
                'Select image file');
            if (path_string~=0)
                
                ResetDisplay(app)

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
                app.ImFInfoArea.Value = printstruct(app.gel_data.imfinfo);

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
                
                app.DeleteBoxButton.Enable = 1;
                app.SaveAnalysisButton.Enable = 1;
                app.OutputButton.Enable = 1;

                temp = load(fullfile(path_string,file_string),'-mat','save_data');
                save_data = temp.save_data;

                % Restore
                app.gel_data = [];
                app.gel_data.image_file_string = save_data.image_file_string;
                app.gel_data.im_data = save_data.im_data;
                app.gel_data.imfinfo = save_data.imfinfo;
                app.ImFInfoArea.Value = printstruct(app.gel_data.imfinfo);

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
            app.SaveAnalysisButton.Enable = 1;
            app.OutputButton.Enable = 1;
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

        % Button pushed function: SaveAnalysisButton
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
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create GelBoxUIFigure and hide until all components are created
            app.GelBoxUIFigure = uifigure('Visible', 'off');
            colormap(app.GelBoxUIFigure, 'parula');
            app.GelBoxUIFigure.Position = [100 100 1530 585];
            app.GelBoxUIFigure.Name = 'GelBox';

            % Create FileControlsPanel
            app.FileControlsPanel = uipanel(app.GelBoxUIFigure);
            app.FileControlsPanel.Title = 'File Controls';
            app.FileControlsPanel.Position = [8 390 143 181];

            % Create LoadImageButton
            app.LoadImageButton = uibutton(app.FileControlsPanel, 'push');
            app.LoadImageButton.ButtonPushedFcn = createCallbackFcn(app, @LoadImageButtonPushed, true);
            app.LoadImageButton.Position = [27 134 90 22];
            app.LoadImageButton.Text = 'Load Image';

            % Create InvertImageButton
            app.InvertImageButton = uibutton(app.FileControlsPanel, 'push');
            app.InvertImageButton.ButtonPushedFcn = createCallbackFcn(app, @InvertImageButtonPushed, true);
            app.InvertImageButton.Position = [27 102 90 22];
            app.InvertImageButton.Text = 'Invert Image';

            % Create LoadAnalysisButton
            app.LoadAnalysisButton = uibutton(app.FileControlsPanel, 'push');
            app.LoadAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @LoadAnalysisButtonPushed, true);
            app.LoadAnalysisButton.Position = [27 71 90 22];
            app.LoadAnalysisButton.Text = 'Load Analysis';

            % Create SaveAnalysisButton
            app.SaveAnalysisButton = uibutton(app.FileControlsPanel, 'push');
            app.SaveAnalysisButton.ButtonPushedFcn = createCallbackFcn(app, @SaveAnalysisButtonPushed, true);
            app.SaveAnalysisButton.Enable = 'off';
            app.SaveAnalysisButton.Position = [27 39 90 23];
            app.SaveAnalysisButton.Text = 'Save Analysis';

            % Create OutputButton
            app.OutputButton = uibutton(app.FileControlsPanel, 'push');
            app.OutputButton.Enable = 'off';
            app.OutputButton.Position = [27 8 90 22];
            app.OutputButton.Text = 'Output';

            % Create AnalysisControlsPanel
            app.AnalysisControlsPanel = uipanel(app.GelBoxUIFigure);
            app.AnalysisControlsPanel.Title = 'Analysis Controls';
            app.AnalysisControlsPanel.Position = [157 391 239 181];

            % Create NumberofBandsDropDownLabel
            app.NumberofBandsDropDownLabel = uilabel(app.AnalysisControlsPanel);
            app.NumberofBandsDropDownLabel.Position = [9 111 99 22];
            app.NumberofBandsDropDownLabel.Text = 'Number of Bands';

            % Create NumberofBandsDropDown
            app.NumberofBandsDropDown = uidropdown(app.AnalysisControlsPanel);
            app.NumberofBandsDropDown.Items = {'1', '2', '3'};
            app.NumberofBandsDropDown.ValueChangedFcn = createCallbackFcn(app, @NumberofBandsDropDownValueChanged, true);
            app.NumberofBandsDropDown.Position = [123 111 101 22];
            app.NumberofBandsDropDown.Value = '1';

            % Create NewBoxButton
            app.NewBoxButton = uibutton(app.AnalysisControlsPanel, 'push');
            app.NewBoxButton.ButtonPushedFcn = createCallbackFcn(app, @NewBoxButtonPushed, true);
            app.NewBoxButton.Position = [9 73 100 22];
            app.NewBoxButton.Text = 'New Box';

            % Create DeleteBoxButton
            app.DeleteBoxButton = uibutton(app.AnalysisControlsPanel, 'push');
            app.DeleteBoxButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteBoxButtonPushed, true);
            app.DeleteBoxButton.Enable = 'off';
            app.DeleteBoxButton.Position = [123 73 101 22];
            app.DeleteBoxButton.Text = 'Delete Box';

            % Create BoxSelectionDropDownLabel
            app.BoxSelectionDropDownLabel = uilabel(app.AnalysisControlsPanel);
            app.BoxSelectionDropDownLabel.Position = [10 33 98 22];
            app.BoxSelectionDropDownLabel.Text = 'Box Selection';

            % Create BoxSelectionDropDown
            app.BoxSelectionDropDown = uidropdown(app.AnalysisControlsPanel);
            app.BoxSelectionDropDown.Items = {};
            app.BoxSelectionDropDown.ValueChangedFcn = createCallbackFcn(app, @BoxSelectionDropDownValueChanged, true);
            app.BoxSelectionDropDown.Placeholder = 'No data';
            app.BoxSelectionDropDown.Position = [124 33 100 22];
            app.BoxSelectionDropDown.Value = {};

            % Create SelectedBoxInformationPanel
            app.SelectedBoxInformationPanel = uipanel(app.GelBoxUIFigure);
            app.SelectedBoxInformationPanel.Title = 'Selected Box Information';
            app.SelectedBoxInformationPanel.Position = [403 390 283 182];

            % Create SelectedBoxInformationArea
            app.SelectedBoxInformationArea = uitextarea(app.SelectedBoxInformationPanel);
            app.SelectedBoxInformationArea.Editable = 'off';
            app.SelectedBoxInformationArea.FontSize = 11;
            app.SelectedBoxInformationArea.Position = [12 11 260 143];

            % Create SelectedBoxOpticalDensitiesPanel
            app.SelectedBoxOpticalDensitiesPanel = uipanel(app.GelBoxUIFigure);
            app.SelectedBoxOpticalDensitiesPanel.Title = 'Selected Box Optical Densities';
            app.SelectedBoxOpticalDensitiesPanel.Position = [698 294 826 278];

            % Create raw_density
            app.raw_density = uiaxes(app.SelectedBoxOpticalDensitiesPanel);
            xlabel(app.raw_density, 'Optical Density')
            ylabel(app.raw_density, 'Pixel')
            app.raw_density.Position = [163 11 254 209];

            % Create background_corrected_raw_density
            app.background_corrected_raw_density = uiaxes(app.SelectedBoxOpticalDensitiesPanel);
            xlabel(app.background_corrected_raw_density, 'Optical Density')
            ylabel(app.background_corrected_raw_density, 'Pixel')
            app.background_corrected_raw_density.YColor = [0.9412 0.9412 0.9412];
            app.background_corrected_raw_density.Position = [401 11 254 209];

            % Create box_inset
            app.box_inset = uiaxes(app.SelectedBoxOpticalDensitiesPanel);
            xlabel(app.box_inset, 'Optical Density')
            ylabel(app.box_inset, 'Pixel')
            app.box_inset.XColor = [0.9412 0.9412 0.9412];
            app.box_inset.YColor = [0.9412 0.9412 0.9412];
            app.box_inset.Position = [4 11 150 209];

            % Create BoxZoomLabel
            app.BoxZoomLabel = uilabel(app.SelectedBoxOpticalDensitiesPanel);
            app.BoxZoomLabel.HorizontalAlignment = 'center';
            app.BoxZoomLabel.WordWrap = 'on';
            app.BoxZoomLabel.Position = [34 224 126 22];
            app.BoxZoomLabel.Text = 'Box Zoom';

            % Create RawOpticalDensityLabel_2
            app.RawOpticalDensityLabel_2 = uilabel(app.SelectedBoxOpticalDensitiesPanel);
            app.RawOpticalDensityLabel_2.HorizontalAlignment = 'center';
            app.RawOpticalDensityLabel_2.WordWrap = 'on';
            app.RawOpticalDensityLabel_2.Position = [238 224 126 22];
            app.RawOpticalDensityLabel_2.Text = 'Raw Optical Density';

            % Create BackgroundCorrectedOpticalDensityLabel_2
            app.BackgroundCorrectedOpticalDensityLabel_2 = uilabel(app.SelectedBoxOpticalDensitiesPanel);
            app.BackgroundCorrectedOpticalDensityLabel_2.HorizontalAlignment = 'center';
            app.BackgroundCorrectedOpticalDensityLabel_2.WordWrap = 'on';
            app.BackgroundCorrectedOpticalDensityLabel_2.Position = [477 221 126 28];
            app.BackgroundCorrectedOpticalDensityLabel_2.Text = 'Background Corrected Optical Density';

            % Create TotalAreaEditFieldLabel
            app.TotalAreaEditFieldLabel = uilabel(app.SelectedBoxOpticalDensitiesPanel);
            app.TotalAreaEditFieldLabel.Position = [661 146 59 22];
            app.TotalAreaEditFieldLabel.Text = 'Total Area';

            % Create TotalAreaField
            app.TotalAreaField = uieditfield(app.SelectedBoxOpticalDensitiesPanel, 'numeric');
            app.TotalAreaField.ValueDisplayFormat = '%.2f';
            app.TotalAreaField.Editable = 'off';
            app.TotalAreaField.HorizontalAlignment = 'center';
            app.TotalAreaField.Position = [732 146 85 22];

            % Create BackgroundAreaLabel
            app.BackgroundAreaLabel = uilabel(app.SelectedBoxOpticalDensitiesPanel);
            app.BackgroundAreaLabel.WordWrap = 'on';
            app.BackgroundAreaLabel.Position = [661 106 81 30];
            app.BackgroundAreaLabel.Text = 'Background Area';

            % Create BackgroundAreaField
            app.BackgroundAreaField = uieditfield(app.SelectedBoxOpticalDensitiesPanel, 'numeric');
            app.BackgroundAreaField.ValueDisplayFormat = '%.2f';
            app.BackgroundAreaField.Editable = 'off';
            app.BackgroundAreaField.HorizontalAlignment = 'center';
            app.BackgroundAreaField.Position = [732 110 85 22];

            % Create GelImageFileInformationPanel
            app.GelImageFileInformationPanel = uipanel(app.GelBoxUIFigure);
            app.GelImageFileInformationPanel.Title = 'Gel Image File Information';
            app.GelImageFileInformationPanel.Position = [8 8 211 370];

            % Create ImFInfoArea
            app.ImFInfoArea = uitextarea(app.GelImageFileInformationPanel);
            app.ImFInfoArea.Editable = 'off';
            app.ImFInfoArea.FontSize = 9;
            app.ImFInfoArea.Position = [8 9 193 334];

            % Create GelImagePanel
            app.GelImagePanel = uipanel(app.GelBoxUIFigure);
            app.GelImagePanel.Title = 'Gel Image';
            app.GelImagePanel.Position = [227 8 459 370];

            % Create gel_image_axes
            app.gel_image_axes = uiaxes(app.GelImagePanel);
            app.gel_image_axes.XTick = [];
            app.gel_image_axes.YTick = [];
            app.gel_image_axes.Box = 'on';
            app.gel_image_axes.Position = [11 14 437 329];

            % Create SelectedBoxFittingPanel
            app.SelectedBoxFittingPanel = uipanel(app.GelBoxUIFigure);
            app.SelectedBoxFittingPanel.Title = 'Selected Box Fitting';
            app.SelectedBoxFittingPanel.Position = [698 8 827 277];

            % Create background_corrected_raw_density_fit
            app.background_corrected_raw_density_fit = uiaxes(app.SelectedBoxFittingPanel);
            xlabel(app.background_corrected_raw_density_fit, 'Optical Density')
            ylabel(app.background_corrected_raw_density_fit, 'Pixel')
            app.background_corrected_raw_density_fit.YColor = [0.9412 0.9412 0.9412];
            app.background_corrected_raw_density_fit.Position = [247 13 254 209];

            % Create raw_density_fit
            app.raw_density_fit = uiaxes(app.SelectedBoxFittingPanel);
            xlabel(app.raw_density_fit, 'Optical Density')
            ylabel(app.raw_density_fit, 'Pixel')
            app.raw_density_fit.Position = [9 13 254 209];

            % Create RawOpticalDensityLabel
            app.RawOpticalDensityLabel = uilabel(app.SelectedBoxFittingPanel);
            app.RawOpticalDensityLabel.HorizontalAlignment = 'center';
            app.RawOpticalDensityLabel.WordWrap = 'on';
            app.RawOpticalDensityLabel.Position = [82 224 126 22];
            app.RawOpticalDensityLabel.Text = 'Raw Optical Density';

            % Create BackgroundCorrectedOpticalDensityLabel
            app.BackgroundCorrectedOpticalDensityLabel = uilabel(app.SelectedBoxFittingPanel);
            app.BackgroundCorrectedOpticalDensityLabel.HorizontalAlignment = 'center';
            app.BackgroundCorrectedOpticalDensityLabel.WordWrap = 'on';
            app.BackgroundCorrectedOpticalDensityLabel.Position = [329 221 126 28];
            app.BackgroundCorrectedOpticalDensityLabel.Text = 'Background Corrected Optical Density';

            % Create DrawFittingCheckBox
            app.DrawFittingCheckBox = uicheckbox(app.SelectedBoxFittingPanel);
            app.DrawFittingCheckBox.Text = 'Draw Fitting';
            app.DrawFittingCheckBox.Position = [513 199 86 22];

            % Create RsquaredLabel
            app.RsquaredLabel = uilabel(app.SelectedBoxFittingPanel);
            app.RsquaredLabel.WordWrap = 'on';
            app.RsquaredLabel.Position = [514 167 72 22];
            app.RsquaredLabel.Text = 'R - squared';

            % Create rsquaredField
            app.rsquaredField = uieditfield(app.SelectedBoxFittingPanel, 'numeric');
            app.rsquaredField.ValueDisplayFormat = '%.2f';
            app.rsquaredField.Editable = 'off';
            app.rsquaredField.HorizontalAlignment = 'center';
            app.rsquaredField.Position = [584 167 69 22];

            % Create BandAreaLabel_1
            app.BandAreaLabel_1 = uilabel(app.SelectedBoxFittingPanel);
            app.BandAreaLabel_1.HorizontalAlignment = 'center';
            app.BandAreaLabel_1.WordWrap = 'on';
            app.BandAreaLabel_1.Position = [513 124 70 28];
            app.BandAreaLabel_1.Text = 'Band Area 1 (Red)';

            % Create BandArea_1
            app.BandArea_1 = uieditfield(app.SelectedBoxFittingPanel, 'numeric');
            app.BandArea_1.ValueDisplayFormat = '%.2f';
            app.BandArea_1.Editable = 'off';
            app.BandArea_1.HorizontalAlignment = 'center';
            app.BandArea_1.Position = [584 128 69 22];

            % Create BandRelativeAreaLabel_1
            app.BandRelativeAreaLabel_1 = uilabel(app.SelectedBoxFittingPanel);
            app.BandRelativeAreaLabel_1.WordWrap = 'on';
            app.BandRelativeAreaLabel_1.Position = [663 117 81 42];
            app.BandRelativeAreaLabel_1.Text = 'Band Relative Area 1 (Red)';

            % Create BandRelativeArea_1
            app.BandRelativeArea_1 = uieditfield(app.SelectedBoxFittingPanel, 'numeric');
            app.BandRelativeArea_1.ValueDisplayFormat = '%.2f';
            app.BandRelativeArea_1.Editable = 'off';
            app.BandRelativeArea_1.HorizontalAlignment = 'center';
            app.BandRelativeArea_1.Position = [749 127 70 22];

            % Create BandAreaLabel_2
            app.BandAreaLabel_2 = uilabel(app.SelectedBoxFittingPanel);
            app.BandAreaLabel_2.HorizontalAlignment = 'center';
            app.BandAreaLabel_2.WordWrap = 'on';
            app.BandAreaLabel_2.Enable = 'off';
            app.BandAreaLabel_2.Position = [511 86 72 28];
            app.BandAreaLabel_2.Text = 'Band Area 2 (Blue)';

            % Create BandArea_2
            app.BandArea_2 = uieditfield(app.SelectedBoxFittingPanel, 'numeric');
            app.BandArea_2.ValueDisplayFormat = '%.2f';
            app.BandArea_2.Editable = 'off';
            app.BandArea_2.HorizontalAlignment = 'center';
            app.BandArea_2.Enable = 'off';
            app.BandArea_2.Position = [584 89 69 22];

            % Create BandRelativeAreaLabel_2
            app.BandRelativeAreaLabel_2 = uilabel(app.SelectedBoxFittingPanel);
            app.BandRelativeAreaLabel_2.WordWrap = 'on';
            app.BandRelativeAreaLabel_2.Enable = 'off';
            app.BandRelativeAreaLabel_2.Position = [663 79 81 42];
            app.BandRelativeAreaLabel_2.Text = 'Band Relative Area 2 (Blue)';

            % Create BandRelativeArea_2
            app.BandRelativeArea_2 = uieditfield(app.SelectedBoxFittingPanel, 'numeric');
            app.BandRelativeArea_2.ValueDisplayFormat = '%.2f';
            app.BandRelativeArea_2.Editable = 'off';
            app.BandRelativeArea_2.HorizontalAlignment = 'center';
            app.BandRelativeArea_2.Enable = 'off';
            app.BandRelativeArea_2.Position = [749 89 70 22];

            % Create BandAreaLabel_3
            app.BandAreaLabel_3 = uilabel(app.SelectedBoxFittingPanel);
            app.BandAreaLabel_3.HorizontalAlignment = 'center';
            app.BandAreaLabel_3.WordWrap = 'on';
            app.BandAreaLabel_3.Enable = 'off';
            app.BandAreaLabel_3.Position = [511 49 72 28];
            app.BandAreaLabel_3.Text = 'Band Area 3 (Green)';

            % Create BandArea_3
            app.BandArea_3 = uieditfield(app.SelectedBoxFittingPanel, 'numeric');
            app.BandArea_3.ValueDisplayFormat = '%.2f';
            app.BandArea_3.Editable = 'off';
            app.BandArea_3.HorizontalAlignment = 'center';
            app.BandArea_3.Enable = 'off';
            app.BandArea_3.Position = [584 53 69 22];

            % Create BandRelativeAreaLabel_3
            app.BandRelativeAreaLabel_3 = uilabel(app.SelectedBoxFittingPanel);
            app.BandRelativeAreaLabel_3.WordWrap = 'on';
            app.BandRelativeAreaLabel_3.Enable = 'off';
            app.BandRelativeAreaLabel_3.Position = [663 42 81 42];
            app.BandRelativeAreaLabel_3.Text = 'Band Relative Area 3 (Green)';

            % Create BandRelativeArea_3
            app.BandRelativeArea_3 = uieditfield(app.SelectedBoxFittingPanel, 'numeric');
            app.BandRelativeArea_3.ValueDisplayFormat = '%.2f';
            app.BandRelativeArea_3.Editable = 'off';
            app.BandRelativeArea_3.HorizontalAlignment = 'center';
            app.BandRelativeArea_3.Enable = 'off';
            app.BandRelativeArea_3.Position = [749 52 70 22];

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

