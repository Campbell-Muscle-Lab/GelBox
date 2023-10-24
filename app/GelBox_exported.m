classdef GelBox_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        GelBoxUIFigure                 matlab.ui.Figure
        FileMenu                       matlab.ui.container.Menu
        LoadImageMenu                  matlab.ui.container.Menu
        LoadAnalysisMenu               matlab.ui.container.Menu
        SaveAnalysisMenu               matlab.ui.container.Menu
        ExportResultsMenu              matlab.ui.container.Menu
        DataAnalysisMenu               matlab.ui.container.Menu
        GelImageFileInformationMenu    matlab.ui.container.Menu
        SelectedBoxInformationMenu     matlab.ui.container.Menu
        SummaryPlotMenu                matlab.ui.container.Menu
        FittingPanel                   matlab.ui.container.Panel
        BandTable                      matlab.ui.control.Table
        NumberofBandsSpinner           matlab.ui.control.Spinner
        NumberofBandsSpinnerLabel      matlab.ui.control.Label
        RsquaredLabel                  matlab.ui.control.Label
        rsquaredField                  matlab.ui.control.NumericEditField
        FittingParametersButton        matlab.ui.control.Button
        DrawFittingCheckBox            matlab.ui.control.CheckBox
        BackgroundCorrectedOpticalDensityLabel  matlab.ui.control.Label
        RawOpticalDensityLabel         matlab.ui.control.Label
        raw_density_fit                matlab.ui.control.UIAxes
        background_corrected_raw_density_fit  matlab.ui.control.UIAxes
        GelImagePanel                  matlab.ui.container.Panel
        BoxSelectionDropDown           matlab.ui.control.DropDown
        BoxSelectionDropDownLabel      matlab.ui.control.Label
        DeleteBoxButton                matlab.ui.control.Button
        NewBoxButton                   matlab.ui.control.Button
        AdjustImageButton              matlab.ui.control.Button
        gel_image_axis                 matlab.ui.control.UIAxes
        OpticalDensitiesPanel          matlab.ui.container.Panel
        MedianFilterSizeSpinner        matlab.ui.control.Spinner
        MedianFilterSizeSpinnerLabel   matlab.ui.control.Label
        ApplyFilterCheckBox            matlab.ui.control.CheckBox
        BackgroundSubtractionDropDown  matlab.ui.control.DropDown
        BackgroundSubtractionDropDownLabel  matlab.ui.control.Label
        DensityValueEditField          matlab.ui.control.NumericEditField
        DensityValueEditFieldLabel     matlab.ui.control.Label
        BackgroundAreaField            matlab.ui.control.NumericEditField
        BackgroundAreaLabel            matlab.ui.control.Label
        BackgroundCorrAreaField        matlab.ui.control.NumericEditField
        BackgroundCorrAreaLabel        matlab.ui.control.Label
        RollingBallSizeSpinner         matlab.ui.control.Spinner
        RollingBallSizeSpinnerLabel    matlab.ui.control.Label
        TotalAreaField                 matlab.ui.control.NumericEditField
        TotalAreaEditFieldLabel        matlab.ui.control.Label
        BackgroundCorrectedOpticalDensityLabel_2  matlab.ui.control.Label
        RawOpticalDensityLabel_2       matlab.ui.control.Label
        BoxZoomLabel                   matlab.ui.control.Label
        box_inset                      matlab.ui.control.UIAxes
        background_corrected_raw_density  matlab.ui.control.UIAxes
        raw_density                    matlab.ui.control.UIAxes
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
        image_adjusted = 0
        single_box_callback
        d
        par_fit_na = 0 % Description
        moving_box = 0 % Description
        background_token = 0 % Description
        filtered_inset = [] % Description
    end

    properties (Access = private)
        ImageFileTextDialog % Description
        SelectedBoxTextDialog % Description
        FittingOptions % Description
        SummaryPlot
        AdjustImage
    end

    methods (Access = public)

        function UpdateDisplay(app)
            if (isfield(app.gel_data.boxes,'box_handle'))
                n = numel(app.gel_data.boxes.box_handle);
                t = size(app.gel_data.fitting.par_update,2);
                if t ~= n
                    for l = t+1:n
                        app.gel_data.fitting.par_update(l) = 0;
                    end
                end
                % Get selected box in control
                control_strings = app.BoxSelectionDropDown.Value;
                selected_box = str2num(control_strings);

                for i=1:n
                    p(i,1:4) = app.gel_data.boxes.box_handle(i).Position;
                end
                w = p(selected_box,3);
                h = p(selected_box,4);
                if ((w~=app.gel_data.old_width)|(h~=app.gel_data.old_height))
                    for i=1:n
                        p(i,3) = w;
                        p(i,4) = h;
                        app.gel_data.boxes.box_handle(i).Position = p(i,:);
                    end
                end

                % Store data in case we need to save it
                for i=1:n
                    app.gel_data.boxes.box_position(i,:) = app.gel_data.boxes.box_handle(i).Position;
                end

                if app.single_box_callback
                    disp_start = selected_box;
                    disp_end = selected_box;
                else
                    disp_start = 1;
                    disp_end = n;
                end

                %                 app.d.box(selected_box).fitting_mode = str2num(app.NumberofBandsDropDown.Value);


                for i = disp_start:disp_end

                    %                     app.d.box(i).fitting_mode = str2num(app.NumberofBandsDropDown.Value);
                    num_of_bands = app.d.box(i).fitting_mode;
                    % Extract position
                    app.d.box(i).position = app.gel_data.boxes.box_handle(i).Position;

                    % Label it
                    set(app.gel_data.boxes.box_label(i),'String',sprintf('%.0f',i));
                    set(app.gel_data.boxes.box_label(i), ...
                        'Position',[app.d.box(i).position(1)+app.d.box(i).position(3)+10 ...
                        app.d.box(i).position(2)-50]);

                    % Calculate profile
                    app.d.box(i).inset = imcrop(app.gel_data.image.im_data, ...
                        app.d.box(i).position);
                    summary_position = [app.d.box(i).position(1)-10, app.d.box(i).position(2)-10, app.d.box(i).position(3)+20, app.d.box(i).position(4)+20];
                    app.d.box(i).summary_inset = imcrop(app.gel_data.image.im_data, ...
                        summary_position);
                    
                    if app.ApplyFilterCheckBox.Value
                        app.MedianFilterSizeSpinner.Enable = 'on';
                        app.MedianFilterSizeSpinnerLabel.Enable = 'on';
                        sz = app.MedianFilterSizeSpinner.Value;
                        app.d.box(i).inset = medfilt2(app.d.box(i).inset,[sz sz],'symmetric');
                        app.gel_data.settings.filtering.median.size(i) = sz;
                        app.filtered_inset(i) = 1;
                    else
                        app.MedianFilterSizeSpinner.Enable = 'off';
                        app.MedianFilterSizeSpinnerLabel.Enable = 'off';
                        app.gel_data.settings.filtering.median.size(i) = 0;
                        app.filtered_inset(i) = 0;
                    end

                    method = app.BackgroundSubtractionDropDown.Value;
%                     m = imcomplement(app.d.box(i).inset);
                    m = (app.d.box(i).inset);

                    x = flipud(mean(imcomplement(m),2));
                    y = 1:size(m,1);
                    y = y';
                    
                    
                    if ~app.background_token(i)
                    switch method
                        case 'Rolling Ball'
                            radius = app.RollingBallSizeSpinner.Value;
                            se = strel('disk', radius, 0);
                            app.gel_data.settings.background.method{i} = method;
                            app.gel_data.settings.background.size(i) = radius;
                            im_close = [];
                            im_close = imclose(app.d.box(i).inset,se);
                            im_close = imcomplement(im_close);
                            mean_imclose = mean(im_close,2);
                            mean_imclose = flipud(mean_imclose);
                            app.gel_data.background(i).x_back = mean_imclose;
                            app.background_token(i) = 1;
                        case 'Linear'
                            app.gel_data.background(i).x_back = linspace(x(1),x(end),numel(x));
                            app.gel_data.background(i).x_back = app.gel_data.background(i).x_back';
                            app.gel_data.settings.background.method{i} = method;
                            app.background_token(i) = 1;
                        case 'Constant Value'
                            cst_val = app.DensityValueEditField.Value;
                            app.gel_data.background(i).x_back = cst_val * ones(numel(x),1);
                            app.gel_data.settings.background.method{i} = method;
                            app.gel_data.settings.background.size(i) = cst_val;
                            app.background_token(i) = 1;
                        case 'Cubic Spline'
                            number_of_elements = numel(x);
                            fraction = 0.1;
                            fraction_number = ceil(fraction*numel(x))
                            
                            ix_1 = 1:fraction_number;
                            ix_2 = number_of_elements:-1:(number_of_elements - fraction_number + 1);
                            
                            c_x = [ix_1 ix_2]
                            c_y = x(c_x)
                            
                            app.gel_data.background(i).x_back = csapi(c_x,c_y,1:numel(x))';
                            
                            app.gel_data.settings.background.method{i} = method;
                            app.background_token(i) = 1;
                        case 'Polyline'
                            
                            diff_x = diff(x);
                            
                            [xout_diff,yout_diff] = intersections(1:numel(diff_x),diff_x,...
                                1:numel(diff_x),zeros(1,numel(diff_x)),1);
                            
                            p1 = ceil(xout_diff(1)) - 1;
                            ix = ceil(xout_diff) - 1;
                            
                            k = 3;
                            [idx,C] = kmeans(xout_diff,k);

                            [~,c_ix] = max(C);
                            
                            p2_cluster = xout_diff(idx == c_ix);
                            p2 = ceil(p2_cluster(1)) - 1;
                            app.gel_data.background(i).x_back = [(x(1:p1-1))'...
                                linspace(x(p1),x(p2),(p2-p1+1))...
                                (x(p2+1:end))']';

                            
                            app.gel_data.settings.background.method{i} = method;
                            app.background_token(i) = 1;
                    end
                    end

                    box_no = str2num(app.BoxSelectionDropDown.Value);
                    if (app.loaded_analysis && (app.new_box <= length(app.gel_data.fitting.par_est))) ...
                            || app.parameters_updated || app.gel_data.fitting.par_update(i) || app.moving_box
                    elseif app.new_box || app.mode_updated || app.par_est_na
                        [par_est,par_con] = EstimateFittingParameters(app,y,x, ...
                            app.gel_data.background(i).x_back,num_of_bands);
                        app.gel_data.fitting.par_est(i) = deal(par_est);
                        app.gel_data.fitting.par_con(i) = deal(par_con);
                    end
                    [x_bands,x_fit,r_squared,par_fit] = ...
                        FitGaussian(app,y,x,app.gel_data.background(i).x_back,i,num_of_bands);
                    
                    app.gel_data.fitting.par_fit(i) = deal(par_fit);
                    fnames = fieldnames(par_fit);

                    app.d.box(i).total_area = simps(y,x);
                    app.d.box(i).background_area = simps(y,app.gel_data.background(i).x_back);
                    app.d.box(i).background_corr_area = simps(y,(x-app.gel_data.background(i).x_back));

                    for j = 1 : num_of_bands
                        app.d.box(i).band_area(j) = simps(y,x_bands(:,j));
                        sum(x_bands(:,j));
                    end

                    % Store data for later

                    app.gel_data.box_data(i) = app.d.box(i);
                    app.gel_data.summary(i).x = x;
                    app.gel_data.summary(i).y = y;
                    app.gel_data.summary(i).x_fit = x_fit;
                    app.gel_data.summary(i).x_back = app.gel_data.background(i).x_back;


                    [~,ix_locs] = sort(app.gel_data.fitting.par_fit(i).peak_location);

                    app.gel_data.summary(i).band = [];
                    app.gel_data.summary(i).area = [];
                    for u = 1 : num_of_bands
                        app.gel_data.summary(i).band(:,u) = x_bands(:,ix_locs(u));
                        app.gel_data.summary(i).area(u) = app.d.box(i).band_area(ix_locs(u));
                    end

                    app.gel_data.summary(i).inset = app.d.box(i).inset;
                    app.gel_data.summary(i).summary_inset = app.d.box(i).summary_inset;
                    app.gel_data.summary(i).r_squared = r_squared;
                    app.gel_data.summary(i).box_position = app.d.box(i).position;

                    % Display
                    if (i==selected_box)
                        app.TotalAreaField.Value = app.d.box(i).total_area;
                        app.BackgroundAreaField.Value = app.d.box(i).background_area;
                        app.BackgroundCorrAreaField.Value = app.d.box(i).background_corr_area;
                        app.rsquaredField.Value = r_squared;
                        center_image_with_preserved_aspect_ratio( ...
                            app.d.box(i).inset, ...
                            app.box_inset,[]);
                        cla(app.raw_density)
                        plot(app.raw_density,x,y,"Color",'k',"LineWidth",2)
                        hold(app.raw_density,"on")
                        plot(app.raw_density,app.gel_data.background(i).x_back,y,'-.m',"LineWidth",2)
                        x_pow = ceil(log10(max(x)));
                        x_tick_rounder = 10^(x_pow - 1);
                        x_t_end = ceil(max(x)/x_tick_rounder)*x_tick_rounder;
                        x_t_mid = round(x_t_end/2);
                        x_ticks = [0 x_t_mid x_t_end];
                        xticks(app.raw_density,x_ticks)
                        app.raw_density.XAxis.Exponent = 0;
                        xlim(app.raw_density,[0 x_t_end]);
                        ylim(app.raw_density,[1 max(y)]);
%                         legend(app.raw_density,'','Baseline', ...
%                             'Location','northeast')


                        xticks(app.raw_density_fit,x_ticks);
                        app.raw_density_fit.XAxis.Exponent = 0;

                        xlim(app.raw_density_fit,[0 x_t_end]);
                        ylim(app.raw_density_fit,[1 max(y)]);
                        % legend(app.raw_density)


                        cla(app.background_corrected_raw_density)
                        plot(app.background_corrected_raw_density, ...
                            x-app.gel_data.background(i).x_back,y,'-.k',"LineWidth",2)
                        hold(app.background_corrected_raw_density,"on")
                        plot(app.background_corrected_raw_density, ...
                            zeros(1,numel(y)),y,'-.m',"LineWidth",2)
%                         legend(app.background_corrected_raw_density,'','Baseline', ...
%                             'Location','northeast')
                        ylim(app.background_corrected_raw_density, ...
                            [1 max(y)]);

                        x_pow = ceil(log10(max(x-app.gel_data.background(i).x_back)));
                        x_tick_rounder = 10^(x_pow - 1);
                        x_t_end = ceil(max(x-app.gel_data.background(i).x_back)/x_tick_rounder)*x_tick_rounder;
                        if min(x-app.gel_data.background(i).x_back') < 0
                            x_tick_rounder = -10^(x_pow - 2)*0.25;
                        else
                            x_tick_rounder = 10^(x_pow - 2)*0.25;
                        end
                        x_t_beginning = ceil(min(x-app.gel_data.background(i).x_back)/x_tick_rounder)*x_tick_rounder;
                        x_t_mid = round(x_t_end/2);
                        app.background_corrected_raw_density.XAxis.Exponent = 0;
                        background_method = app.BackgroundSubtractionDropDown.Value;
                        if strcmp(background_method,'Rolling Ball')
                            xlim(app.background_corrected_raw_density,[0 x_t_end])
                            xlim(app.background_corrected_raw_density_fit,[0 x_t_end])
                            x_ticks = [0 x_t_mid x_t_end];

                        else
                            xlim(app.background_corrected_raw_density,[x_t_beginning x_t_end])
                            xlim(app.background_corrected_raw_density_fit,[x_t_beginning x_t_end])
                            x_ticks = [x_t_beginning x_t_mid x_t_end];

                        end
                        xticks(app.background_corrected_raw_density,x_ticks)
                        ylim(app.background_corrected_raw_density,[1 max(y)]);
                        xticks(app.background_corrected_raw_density_fit,x_ticks)
                        app.background_corrected_raw_density_fit.XAxis.Exponent = 0;
                        ylim(app.background_corrected_raw_density_fit,[1 max(y)]);

                        color = parula(num_of_bands);

                        cla(app.raw_density_fit)
                        cla(app.background_corrected_raw_density_fit)



                        for j = 1 : num_of_bands
                            patch(app.background_corrected_raw_density_fit, ...
                                x_bands(:,ix_locs(j)), ...
                                y,color(j,:),'FaceAlpha',0.65, ...
                                'EdgeColor',color(j,:),'EdgeAlpha',0.25, ...
                                'LineWidth',2)
                        end

                        f_color = '#1aff00';

                        hold(app.background_corrected_raw_density_fit,"on")
                        plot(app.background_corrected_raw_density_fit, ...
                            x-app.gel_data.background(i).x_back,y,'-.k',"LineWidth",2)
                        plot(app.background_corrected_raw_density_fit, ...
                            x_fit,y,':',"LineWidth",2,"Color",f_color)

                        hold(app.raw_density_fit,"on")

                        plot(app.raw_density_fit, ...
                            x,y,"Color",'k',"LineWidth",2)
                        plot(app.raw_density_fit, ...
                            app.gel_data.background(i).x_back+x_fit,y,':',"LineWidth",2,"Color",f_color)

                        ylim(app.raw_density_fit, ...
                            [1 max(y)]);
                        ylim(app.background_corrected_raw_density_fit, ...
                            [1 max(y)]);

                        for count = 1 : num_of_bands
                            bt.band_no(count,:) = count;
                            bt.color{count,:} = '';
                            bt.area(count,:) = app.gel_data.summary(i).area(count);
                            bt.relative_area(count,:) = app.gel_data.summary(i).area(count)/sum(app.gel_data.summary(i).area);
                        end

                        t = struct2table(bt);
                        app.BandTable.Data = t;

                        for count = 1 : num_of_bands
                            s = uistyle("BackgroundColor",color(count,:));
                            addStyle(app.BandTable,s,"cell",[count 2])
                        end
                    end
                end
                app.gel_data.d_box = app.d;
            end
        end

        function UpdateFittingOptions(app)
            app.parameters_updated = 1;
            selected_box = str2num(app.BoxSelectionDropDown.Value);
            app.gel_data.fitting.par_update(selected_box) = 1;
            app.single_box_callback = 1;
            UpdateDisplay(app)
            app.parameters_updated = 0;
            app.single_box_callback = 0;

        end

        function [par_est,par_con] = EstimateFittingParameters(app,x,y,y_back, ...
                no_of_bands)

            peaks = find_peaks('x',x, ...
                'y',y, ...
                'min_rel_delta_y',0.05, ...
                'min_x_index_spacing',2);

            % The order of the parameters is as follows:
            %       1) peak_location
            %       2) amplitude
            %       3) width_parameter
            %       4) skew_parameter

            par_est.band_no = (1 : no_of_bands)';
            par_con.band_no = (1 : no_of_bands)';

            if numel(peaks.max_indices) == no_of_bands
                for i = 1 : no_of_bands
                    par_est.peak_location(i,:) = peaks.max_indices(i);
                end
            else
                for i = 1 : no_of_bands
                    par_est.peak_location(i,:) = ...
                        (i/(1+no_of_bands))*length(x);
                end
            end

            target = y;
            target = target - y_back;

            [max_value,~] = max(target);

            half_distance = 0.1*length(x);
            alfa_estimate = -log(0.5)/(half_distance^2);

            for i = 1 : no_of_bands
                par_est.amplitude(i,:) = max_value;
            end

            % The first width and skewness parameters are assigned,
            % and the following are defined as an offset.

            par_est.width_parameter(1,:) = alfa_estimate;
            par_est.skew_parameter(1,:) = 1;

            for i = 2 : no_of_bands
                par_est.width_parameter(i,:) = 0;
                par_est.skew_parameter(i,:) = 0;
            end

            for i = 1 : no_of_bands
                par_con.peak_location(i,:) = false;
                par_con.amplitude(i,:) = false;
                if i > 1
                    par_con.width_parameter(i,:) = true;
                    par_con.skew_parameter(i,:) = true;
                else
                    par_con.width_parameter(i,:) = false;
                    par_con.skew_parameter(i,:) = false;
                end
            end
        end

        function ImageAdjusted(app)
            if ~isempty(app.gel_data.image.adjusted_image)
                app.gel_data.image.im_data = app.gel_data.image.adjusted_image;
                app.image_adjusted = 1;
                center_image_with_preserved_aspect_ratio( ...
                    app.gel_data.image.im_data, ...
                    app.gel_image_axis,[]);
            end
            UpdateDisplay(app)
        end

        function [y_bands, y_fit,r_squared,par_fit] = FitGaussian(app, x, y, y_back, box_no,no_of_bands)
            par_est = struct();
            par_con = struct();

            box_pars = struct();
            box_pars = app.gel_data.fitting.par_est(box_no);

            box_cons = struct();
            box_cons = app.gel_data.fitting.par_con(box_no);

            names = fieldnames(box_pars);
            for m = 1 : numel(names)
                par_est.(names{m}) = box_pars.(names{m});
                par_con.(names{m}) = box_cons.(names{m});
            end

            m = 1;

            no_of_parameters = (numel(names)-1) * numel(par_est.band_no);

            A_constraints = zeros(1,no_of_parameters);
            A_constraints(1) = 1;
            B_constants = [0];

            lower_bounds=zeros(1,no_of_parameters);
            upper_bounds=Inf*ones(1,no_of_parameters);
            var_per_band = 4;

            upper_bounds(1:var_per_band:end) = numel(y);

            % Unpack the parameter estimate structure;
            for i = 1 : no_of_bands
                for j = 2 : numel(names)
                    par(m) = par_est.(names{j})(i);
                    if par_con.(names{j})(i)
                        lower_bounds(m) = par_est.(names{j})(i);
                        upper_bounds(m) = par_est.(names{j})(i);
                    end
                    m = m + 1;
                end
            end

            target = y;

            target = target - y_back;

            j = 1;
            e = [];

            if no_of_parameters > var_per_band
                for p = 2 : no_of_bands
                    if ~par_con.width_parameter(p)
                        lower_bounds((p*var_per_band)-1) = -inf;
                    end
                    if ~par_con.skew_parameter(p)
                        lower_bounds(p*var_per_band) = -inf;
                    end
                end
            end

            opts=optimset('fminsearch');
            opts.Display='off';
            opts.MaxIter=1000;
            opts.MaxFunEvals=10000;

            [p_result,fval,exitflag,output]= ...
                fminsearchbnd(@profile_error_ngaussian,par, ...
                lower_bounds,upper_bounds,opts);
            r_squared = calculate_r_squared(y,y_fit+y_back);
            error = (y - (y_fit+y_back)).^2;
            error = sum(error);

            par_fit.band_no = [1:no_of_bands]';

            count = 1;
            for n = 1 : no_of_bands
                for m = 2 : numel(names)
                    par_fit.(names{m})(n) = p_result(count);
                    count = count + 1;
                end
            end

            function trial_e = profile_error_ngaussian(par)

                [y_bands, y_fit] = calculate_nprofile(x,par);

                no_of_par = 4;
                no_of_bands = numel(par)/no_of_par;

                e(j) = 0;
                for i  = 1 : numel(target)
                    e(j) = e(j) + (y_fit(i) - target(i))^2;
                end

                trial_e = e(end);

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
            function [y_bands,y_fit] = calculate_nprofile(x,par)

                no_of_par = 4;
                no_of_bands = numel(par)/no_of_par;

                par_map = 1:no_of_par;
                par_map = repmat(par_map,no_of_bands,1);

                for k = 2 : no_of_bands
                    par_map(k,:) = par_map(k,:) + (k-1)*no_of_par;
                end

                for u = 1 : no_of_bands
                    if u == 1
                        y_bands(:,u) = skewed_Gaussian(x, ...
                            par(par_map(u,1)),par(par_map(u,2)), ...
                            par(par_map(u,3)),par(par_map(u,4)));
                    else
                        y_bands(:,u) = skewed_Gaussian(x, ...
                            par(par_map(u,1)),par(par_map(u,2)), ...
                            par(par_map(1,3)) + par(par_map(u,3)), ...
                            par(par_map(1,4)) + par(par_map(u,4)));
                    end
                end

                y_fit = 0;
                for t = 1 : no_of_bands
                    y_fit = y_fit + y_bands(:,t);
                end

            end
            function y=skewed_Gaussian(x,x0,A,gamma,skew1)
                offset = zeros(length(x),1);
                offset((x-x0)>0) = skew1*(x((x-x0)>0)-x0);
                y=  A*exp(-gamma*(((x-x0)+offset).^2));
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
            app.BandTable.Data = [];

            % Reset number of bands
            app.NumberofBandsSpinner.Value = 1;

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

            addpath(genpath('utilities'));
            movegui(app.GelBoxUIFigure,'center')
            app.gel_data.boxes = [];
            
            app.MedianFilterSizeSpinner.Enable = 'off';
            app.MedianFilterSizeSpinnerLabel.Enable = 'off';


            colormap(app.GelBoxUIFigure, 'gray');
            app.fitting_options.shared_shape = true(1,1);
            app.fitting_options.shared_skewness = true(1,1);
            disableDefaultInteractivity(app.gel_image_axis)

        end

        % Menu selected function: LoadImageMenu
        function LoadImageButtonPushed(app, event)
            [file_string,path_string]=uigetfile2( ...
                {'*.tif','TIF';'*.png','PNG';'*.bmp','BMP';'*.gif','GIF';'*.jpeg','JPEG';'*.jpeg2000','JPEG2000';'*.pbm','PBM';'*.pgm','PGM'}, ...
                'Select Image File');
            if (path_string~=0)

                ResetDisplay(app)
                app.DataAnalysisMenu.Enable = 1;
                app.gel_data = [];
                app.gel_data.fitting.par_update = 0;
                app.gel_data.boxes = [];
                app.d = [];
                app.filtered_inset = [];
                app.background_token = [];

                app.gel_data.invert_status = 0;
                app.gel_data.image.image_file_string = fullfile(path_string,file_string);
                app.gel_data.image.im_data = imread(app.gel_data.image.image_file_string);
                if (ndims(app.gel_data.image.im_data)==3)
                    app.gel_data.image.im_data = rgb2gray(app.gel_data.image.im_data);
                end
                app.gel_data.image.original_image = app.gel_data.image.im_data;
                center_image_with_preserved_aspect_ratio( ...
                    app.gel_data.image.im_data, ...
                    app.gel_image_axis,[]);

                app.gel_data.image.imfinfo = imfinfo(app.gel_data.image.image_file_string);

            end

        end

        % Menu selected function: LoadAnalysisMenu
        function LoadAnalysisButtonPushed(app, event)

            [file_string,path_string] = uigetfile2( ...
                {'*.gbx','GelBox file';'*.gdf','Gel data file'},'Select GBX File To Load Analysis');

            if (path_string~=0)
                % Delete any old boxes
                if (isfield(app.gel_data.boxes,'box_handle'))
                    n = numel(app.gel_data.boxes.box_handle);
                    for i=1:n
                        delete(app.gel_data.boxes.box_handle(i));
                    end
                end
                app.DeleteBoxButton.Enable = 1;
                app.DataAnalysisMenu.Enable = 1;

                temp = load(fullfile(path_string,file_string),'-mat','save_data');
                save_data = temp.save_data;
                
                if ~isfield(save_data,'settings')
                    new_save_data.boxes.box_position = save_data.box_position;
                    new_save_data.fitting.fitting_mode = save_data.fitting_mode;
                    new_save_data.fitting.par_est = save_data.par_est;
                    new_save_data.fitting.par_fit = save_data.par_fit;
                    new_save_data.fitting.par_con = save_data.par_con;
                    new_save_data.image.image_file_string = save_data.image_file_string;
                    new_save_data.image.im_data = save_data.im_data;
                    new_save_data.image.imfinfo = save_data.imfinfo;
                    new_save_data.image.original_image = save_data.original_image;
                    new_save_data.image.adjusted_image = save_data.adjusted_image;
                    new_save_data.settings.image_adjustments = save_data.image_adjustments;
                    new_save_data.settings.background.method = {};
                    new_save_data.settings.background.size = [];
                    for i = 1 : numel(new_save_data.fitting.par_fit)
                        new_save_data.settings.background.method{i} = save_data.background_method;
                        new_save_data.settings.background.size(i) = save_data.background_size;
                    end
                    save_data = [];
                    save_data = new_save_data;
                end

                % Restore
                app.gel_data = [];
                app.gel_data.boxes = [];
                app.filtered_inset = [];
                app.background_token = [];
                app.d=[];
                save_fields = {'image','fitting','settings'};

                for i = 1 : numel(save_fields)
                    app.gel_data.(save_fields{i}) = save_data.(save_fields{i});
                end

                app.new_box = 0;


                center_image_with_preserved_aspect_ratio( ...
                    app.gel_data.image.im_data, ...
                    app.gel_image_axis,[]);

                n=size(save_data.boxes.box_position,1);
                                
                if ~isfield(app.gel_data.settings,'filtering')
                    app.gel_data.settings.filtering.median.size = 3*ones(1,n);
                end
                control_strings = [];
                app.gel_data.fitting.par_update = zeros(1,n);

                for i=1:n
                    app.gel_data.boxes.box_handle(i) = images.roi.Rectangle(app.gel_image_axis, ...
                        'Position',save_data.boxes.box_position(i,:));
                    control_strings{i} = sprintf('%.0f',i);
                    app.filtered_inset(i) = 0;
                    app.background_token(i) = 0;
                    if isfield(save_data.fitting,'fitting_mode')
                        app.d.box(i).fitting_mode = save_data.fitting.fitting_mode(i);
                    else
                        app.d.box(i).fitting_mode = band_no;
                    end
                end

                app.BoxSelectionDropDown.Items = control_strings;
                app.BoxSelectionDropDown.Value = control_strings{1};

                for i=1:n
                    app.gel_data.boxes.box_handle(i).FaceAlpha = 0;
                    if (i~=1)
                        app.gel_data.boxes.box_handle(i).Color = [1 0 0];
                        app.gel_data.boxes.box_handle(i).InteractionsAllowed = 'none';
                    else
                        app.gel_data.boxes.box_handle(i).Color = [0 1 0];
                        app.gel_data.boxes.box_handle(i).InteractionsAllowed = 'all';

                    end

                    p = app.gel_data.boxes.box_handle(i).Position;
                    app.gel_data.boxes.box_label(i) = text(p(1)+p(3),p(2)-50,sprintf('%.0f',i), ...
                        'Parent',app.gel_image_axis,'FontWeight',"bold","FontSize",18);

                    app.gel_data.old_width = p(3);
                    app.gel_data.old_height = p(4);

                    i=i;
                    addlistener(app.gel_data.boxes.box_handle(i),"MovingROI",@(src,evt) new_box_position2(evt));
                end
                app.BackgroundSubtractionDropDown.Value = save_data.settings.background.method(1);
                app.RollingBallSizeSpinner.Value = save_data.settings.background.size(1);
                app.NumberofBandsSpinner.Value = save_data.fitting.fitting_mode(1);
                if app.gel_data.settings.filtering.median.size(1)
                    app.ApplyFilterCheckBox.Value = 1;
                    app.MedianFilterSizeSpinner.Enable = ' on';
                    app.MedianFilterSizeSpinnerLabel.Enable = 'on';
                    app.MedianFilterSizeSpinner.Value = app.gel_data.settings.filtering.median.size(1);
                else
                    app.ApplyFilterCheckBox.Value = 0;
                    app.MedianFilterSizeSpinner.Enable = ' off';
                    app.MedianFilterSizeSpinnerLabel.Enable = 'off';    
                end
                app.gel_data.settings = save_data.settings;
                drawnow;
                UpdateDisplay(app)
                app.par_est_na = 0;
            end



            % Nested function
            function new_box_position2(evt);
                if (isfield(app.gel_data.boxes,'box_position'))
                    box_position = app.gel_data.boxes.box_position;
                    [r,c]=size(box_position);
                    if (r>=n)&(~isequal(box_position(n,:),evt.CurrentPosition))
                        app.new_box = n;
                        old_size = box_position(n,3:4);
                        current_size = evt.CurrentPosition(3:4);
                        if isequal(old_size,current_size)
                            app.single_box_callback = 1;
                            app.background_token(n) = 0;
                            app.filtered_inset(n) = 0;
                        else
                            app.background_token(1:numel(app.gel_data.boxes.box_handle)) = 0;
                            app.filtered_inset(1:numel(app.gel_data.boxes.box_handle)) = 0;
                        end
                        app.moving_box = 1;
                        UpdateDisplay(app);
                        app.single_box_callback = 0;
                        app.moving_box = 0;

                    end
                else
                    app.new_box = n;
                    app.single_box_callback = 1;
                    app.moving_box = 1;
                    app.background_token(n) = 0;
                    app.filtered_inset(n) = 0;
                    UpdateDisplay(app);
                    app.single_box_callback = 0;
                    app.moving_box = 0;

                end
            end
        end

        % Button pushed function: NewBoxButton
        function NewBoxButtonPushed(app, event)
            if (~isfield(app.gel_data.boxes,'box_handle'))
                app.DeleteBoxButton.Enable = 1;
                n=1;
                app.gel_data.boxes.box_handle(n) = drawrectangle(app.gel_image_axis);
                p = app.gel_data.boxes.box_handle(n).Position;
                app.gel_data.old_width = p(3);
                app.gel_data.old_height = p(4);
            else
                n = 1 + numel(app.gel_data.boxes.box_handle);
                p = app.gel_data.boxes.box_handle(n-1).Position;

                app.gel_data.boxes.box_handle(n) = images.roi.Rectangle(app.gel_image_axis, ...
                    'Position',p + [150,0,0,0]);
                for i=1:(n-1)
                    app.gel_data.boxes.box_handle(i).InteractionsAllowed = 'none';
                end

            end

            addlistener(app.gel_data.boxes.box_handle(n),"MovingROI", ...
                @(src,evt) new_box_position(evt));

            % Set color to last box
            app.gel_data.boxes.box_handle(n).Color = [0 1 0];
            app.gel_data.boxes.box_handle(n).FaceAlpha = 0;
            for i=1:(n-1)
                app.gel_data.boxes.box_handle(i).Color = [1 0 0];
            end

            % Add in a label
            p = app.gel_data.boxes.box_handle(n).Position;
            app.gel_data.boxes.box_label(n) = text(app.gel_image_axis, ...
                p(1)+p(3),p(2)-50,sprintf('%.0f',n),'FontWeight',"bold","FontSize",18);

            % Update zoom control
            for i=1:n
                control_strings{i}=sprintf('%.0f',i);
            end

            app.BoxSelectionDropDown.Items = control_strings;
            app.BoxSelectionDropDown.Value = control_strings{n};

            app.DeleteBoxButton.Enable = 1;
            app.SelectedBoxInformationMenu.Enable = 1;
            app.new_box = n;
            app.single_box_callback = 1;
            app.d.box(n).fitting_mode = app.NumberofBandsSpinner.Value;
            app.filtered_inset(n) = 0;
            app.background_token(n) = 0;
            app.ApplyFilterCheckBox.Value = 0;
            app.MedianFilterSizeSpinner.Value = 3;
            app.MedianFilterSizeSpinner.Enable = 'off';
            app.MedianFilterSizeSpinnerLabel.Enable = 'off';
            UpdateDisplay(app)
            app.single_box_callback = 0;

            function new_box_position(evt);
                if (isfield(app.gel_data.boxes,'box_position'))
                    box_position = app.gel_data.boxes.box_position;
                    [r,c]=size(box_position);
                    if (r>=n)&(~isequal(box_position(n,:),evt.CurrentPosition))
                        app.new_box = n;
                        old_size = box_position(n,3:4);
                        current_size = evt.CurrentPosition(3:4);
                        if isequal(old_size,current_size)
                            app.single_box_callback = 1;
                            app.background_token(n) = 0;
                            app.filtered_inset(n) = 0;
                        else
                            app.background_token(1:numel(app.gel_data.boxes.box_handle)) = 0;
                            app.filtered_inset(1:numel(app.gel_data.boxes.box_handle)) = 0;
                        end
                        app.moving_box = 1;
                        UpdateDisplay(app);
                        app.single_box_callback = 0;
                        app.moving_box = 0;

                    end
                else
                    app.new_box = n;
                    app.single_box_callback = 1;
                    app.moving_box = 1;
                    app.background_token(n) = 0;
                    app.filtered_inset(n) = 0;
                    UpdateDisplay(app);
                    app.single_box_callback = 0;
                    app.moving_box = 0;

                end
            end

        end

        % Value changed function: NumberofBandsSpinner
        function NumberofBandsSpinnerValueChanged(app, event)
            value = app.NumberofBandsSpinner.Value;
            app.mode_updated = 1;
            for u = 1:numel(app.gel_data.fitting.par_update)
                app.gel_data.fitting.par_update(u) = 0;
            end

            selected_box = str2num(app.BoxSelectionDropDown.Value);
            app.d.box(selected_box).fitting_mode = app.NumberofBandsSpinner.Value;
            app.single_box_callback = 1;
            if app.loaded_analysis
                app.loaded_analysis = 0;
            end
            UpdateDisplay(app)
            app.single_box_callback = 0;
            app.mode_updated = 0;

        end

        % Menu selected function: SaveAnalysisMenu
        function SaveAnalysisButtonPushed(app, event)
            save_data.image.image_file_string = app.gel_data.image.image_file_string;
            try
                save_data.boxes.box_position = app.gel_data.boxes.box_position;
            catch
                warndlg('The analysis boxes are not available.')
                return
            end
            save_fields = {'image','fitting','settings'};
            
            for i = 1 : numel(save_fields)
            save_data.(save_fields{i}) = app.gel_data.(save_fields{i});
            end

            number_of_boxes = size(save_data.boxes.box_position,1);

            for i = 1:number_of_boxes
                save_data.fitting.fitting_mode(i) = app.d.box(i).fitting_mode;
            end
            
            [file_string,path_string] = uiputfile2( ...
                {'*.gbx','GelBox file';'*.gdf','Gel data file'},'Select File Name To Save Analysis');

            if (path_string~=0)
                save(fullfile(path_string,file_string),'save_data');
                json_file_name = strrep(fullfile(path_string,file_string),'.gbx','.json');
                savejson('GelBox Settings',app.gel_data.settings,json_file_name);
            end
        end

        % Value changed function: BoxSelectionDropDown
        function BoxSelectionDropDownValueChanged(app, event)
            selected_box = str2num(app.BoxSelectionDropDown.Value);
            control_strings = app.BoxSelectionDropDown.Items;

            n = numel(control_strings);
            for i=1:n
                if (i~=selected_box)
                    app.gel_data.boxes.box_handle(i).Color = [1 0 0];
                    app.gel_data.boxes.box_handle(i).InteractionsAllowed = 'none';
                else
                    app.gel_data.boxes.box_handle(i).Color = [0 1 0];
                    app.gel_data.boxes.box_handle(i).InteractionsAllowed = 'all';
                end
            end

            value = app.d.box(selected_box).fitting_mode;
            app.NumberofBandsSpinner.Value = value;
            if app.gel_data.settings.filtering.median.size(selected_box)
                app.ApplyFilterCheckBox.Value = 1;
                app.MedianFilterSizeSpinner.Enable = 'on';
                app.MedianFilterSizeSpinnerLabel.Enable = 'on';
                app.MedianFilterSizeSpinner.Value = app.gel_data.settings.filtering.median.size(selected_box);
            else
                app.ApplyFilterCheckBox.Value = 0;
                app.MedianFilterSizeSpinner.Enable = 'off';
                app.MedianFilterSizeSpinnerLabel.Enable = 'off';
            end
            app.BackgroundSubtractionDropDown.Value = app.gel_data.settings.background.method(selected_box);
            background_method = app.BackgroundSubtractionDropDown.Value;
            
            switch background_method
                case 'Rolling Ball'
                    app.RollingBallSizeSpinner.Visible = 1;
                    app.RollingBallSizeSpinnerLabel.Visible = 1;
                    app.DensityValueEditField.Visible = 0;
                    app.DensityValueEditFieldLabel.Visible = 0;
                    app.RollingBallSizeSpinner.Value = app.gel_data.settings.background.size(selected_box);
                case 'Linear'
                    app.RollingBallSizeSpinner.Visible = 0;
                    app.RollingBallSizeSpinnerLabel.Visible = 0;
                    app.DensityValueEditField.Visible = 0;
                    app.DensityValueEditFieldLabel.Visible = 0;
                case 'Constant Value'
                    app.DensityValueEditField.Visible = 1;
                    app.DensityValueEditFieldLabel.Visible = 1;
                    app.RollingBallSizeSpinner.Visible = 0;
                    app.RollingBallSizeSpinnerLabel.Visible = 0;
                    app.DensityValueEditField.Value = app.gel_data.settings.background.size(selected_box);
            end
            app.single_box_callback = 1;
            UpdateDisplay(app)
            app.single_box_callback = 0;
        end

        % Button pushed function: DeleteBoxButton
        function DeleteBoxButtonPushed(app, event)
            selected_box = str2num(app.BoxSelectionDropDown.Value);
            control_strings = app.BoxSelectionDropDown.Items;
            if numel(app.gel_data.boxes) == 1
                delete(app.gel_data.boxes.box_handle(selected_box))
                rm = {'box_handle','box_label','box_position'};
                app.gel_data.boxes = rmfield(app.gel_data,'box_handle');
                app.BoxSelectionDropDown.Items = {};
                app.BoxSelectionDropDown.Value = {};
                app.BoxSelectionDropDown.Placeholder = 'No Data';
            else
            delete(app.gel_data.boxes.box_handle(selected_box))
            delete(app.gel_data.boxes.box_label(selected_box))
            field_names = fieldnames(app.gel_data.fitting.par_est);

            for u = 1 : numel(field_names)
                app.gel_data.fitting.par_est(selected_box).(field_names{u}) = [];
            end

            app.gel_data.boxes.box_handle(selected_box) = [];
            app.gel_data.boxes.box_label(selected_box) = [];
            app.gel_data.boxes.box_position(selected_box) = [];

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

            new_selected_box = str2num(app.BoxSelectionDropDown.Value);

            for i=1:n
                delete(app.gel_data.boxes.box_label(i));
                if (i~=new_selected_box)
                    app.gel_data.boxes.box_handle(i).Color = [1 0 0];
                    app.gel_data.boxes.box_handle(i).InteractionsAllowed = 'none';
                else
                    app.gel_data.boxes.box_handle(i).Color = [0 1 0];
                    app.gel_data.boxes.box_handle(i).InteractionsAllowed = 'all';
                end
            end

            app.gel_data.boxes.box_label = [];
            for i = 1:n
                p = app.gel_data.boxes.box_handle(i).Position;
                app.gel_data.boxes.box_label(i) = text(p(1)+p(3),p(2)-50,sprintf('%.0f',i), ...
                    'Parent',app.gel_image_axis,'FontWeight',"bold",'FontSize',18);

                app.gel_data.old_width = p(3);
                app.gel_data.old_height = p(4);
            end

            app.single_box_callback = 1;
            UpdateDisplay(app)
            app.single_box_callback = 0;
            end
        end

        % Menu selected function: ExportResultsMenu
        function OutputButtonPushed(app, event)
            o = [];
            n = numel(app.gel_data.boxes.box_handle);
            for i = 1 : n
                fit_mode(i) = app.gel_data.box_data(i).fitting_mode;
            end
            max_num_bands = max(fit_mode);
            for i=1:n
                o.box(i) = i;
                o.image_file{i} = app.gel_data.image.image_file_string;
                o.total_area(i) = app.gel_data.box_data(i).total_area;
                o.background_method{i} = app.gel_data.settings.background.method{i};
                if strcmp(o.background_method{i},'Rolling Ball')
                    o.background_size(i) = app.gel_data.settings.background.size(i);
                elseif strcmp(o.background_method{i}, 'Constant Value')
                    o.background_size(i) = app.gel_data.settings.background.size(i);
                end
                o.median_filter_size(i) = app.gel_data.settings.filtering.median.size(i);
                o.background_area(i) = app.gel_data.box_data(i).background_area;
                num_of_bands = numel(app.gel_data.box_data(i).band_area);

                o.band_left(i) = app.gel_data.box_data(i).position(1);
                o.band_top(i) = app.gel_data.box_data(i).position(2);
                o.band_width(i) = app.gel_data.box_data(i).position(3);
                o.band_height(i) = app.gel_data.box_data(i).position(4);
                o.fitting_mode{i} = app.gel_data.box_data(i).fitting_mode;
                o.num_of_bands(i) = num_of_bands;
                o.r_squared(i) = app.gel_data.summary(i).r_squared;

                for c = 1 : num_of_bands

                    name = sprintf('band_area_%i',c);
                    o.(name)(i) = app.gel_data.summary(i).area(c);
                    
                    
%                     l = length(o.(name)(c));
                end
                
                if num_of_bands < max_num_bands
                    for count = num_of_bands + 1 : max_num_bands
                        name = sprintf('band_area_%i',count);
                        o.(name)(i) = 0;
                    end
                end
                
                
            end

            [file_string,path_string] = uiputfile2( ...
                {'*.xlsx','Excel file'},'Enter Excel File Name For Analysis Results');
            layout_table = [];
            if (path_string~=0)
            d_out = o;
            names = fieldnames(d_out);
            [file_string_layout,path_string_layout] = uigetfile2( ...
                {'*.xlsx','Excel file'},'Select The Gel Layout File For The Results Excel File');
            
            if (path_string_layout~=0)
                
                layout_file = fullfile(path_string_layout,file_string_layout);
                layout_table = readtable(layout_file);
                layout_names = layout_table.Properties.VariableNames;
                
                while ~any(strcmpi(names,'box'))
                    message = ["Gel layout does not have ""Box"" as one of the columns.","Please make sure to add ""Box"" to Excel file."];
                    selection = uiconfirm(app.GelBoxUIFigure,message,"Reload Gel Layout",...
                        "Options",["Reload Gel Layout","Skip"],"DefaultOption",1,"CancelOption",2,'Icon','error');
                    
                    if ~strcmp(selection,'Skip')
                        [file_string_layout,path_string_layout] = uigetfile2( ...
                            {'*.xlsx','Excel file'},'Select The Gel Layout File For The Results Excel File');
                        
                        layout_file = fullfile(path_string_layout,file_string_layout);
                        layout_table = readtable(layout_file);
                    else
                        
                        layout_table = [];
                        return
                    end
                    
                end
            end
            
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
            
            if ~isempty(layout_table)
                for i = 1 : numel(layout_names)
                    if strcmpi(layout_names{i},'box')
                        layout_table = renamevars(layout_table,layout_names{i},'box');
                    end
                end
                T = innerjoin(layout_table,struct2table(d_out));
                writetable(T,output_file,'Sheet','Summary')
            else
                writetable(struct2table(d_out),output_file,'Sheet','Summary')
            end

            summary = app.gel_data.summary;

            for i = 1:n
                i_name = sprintf('box_%i',i);
            end
            remf = {'area','r_squared','inset','summary_inset','box_position'};
            summary = rmfield(summary,remf);

            for box = 1 : size(summary,2)
                for i = 1 : max(o.num_of_bands)
                    name = sprintf('band_%i',i);
                    if size(summary(box).band,2) < i
                        summary(box).(name) = NaN * ones(1,size(summary(box).band(:,1),1));
                    else
                        summary(box).(name) = summary(box).band(:,i);
                    end
                end
            end

            summary = rmfield(summary,'band');
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

        % Value changed function: DrawFittingCheckBox
        function DrawFittingCheckBoxValueChanged(app, event)
            value = app.DrawFittingCheckBox.Value;
            if value
                app.single_box_callback = 1;
                UpdateDisplay(app)
                app.single_box_callback = 0;
            end
        end

        % Button pushed function: AdjustImageButton
        function AdjustImageButtonPushed(app, event)
            app.AdjustImage = AdjustImageWindow(app);
        end

        % Value changed function: RollingBallSizeSpinner
        function RollingBallSizeSpinnerValueChanged(app, event)
            value = app.RollingBallSizeSpinner.Value;
            box = str2num(app.BoxSelectionDropDown.Value);
            app.background_token(box) = 0;
            app.single_box_callback = 1;
            UpdateDisplay(app)
            app.single_box_callback = 0;
        end

        % Value changed function: BackgroundSubtractionDropDown
        function BackgroundSubtractionDropDownValueChanged(app, event)
            value = app.BackgroundSubtractionDropDown.Value;
            switch value
                case 'Rolling Ball'
                    app.RollingBallSizeSpinner.Visible = 1;
                    app.RollingBallSizeSpinnerLabel.Visible = 1;
                    app.DensityValueEditField.Visible = 0;
                    app.DensityValueEditFieldLabel.Visible = 0;
                case 'Linear'
                    app.RollingBallSizeSpinner.Visible = 0;
                    app.RollingBallSizeSpinnerLabel.Visible = 0;
                    app.DensityValueEditField.Visible = 0;
                    app.DensityValueEditFieldLabel.Visible = 0;
                case 'Constant Value'
                    app.DensityValueEditField.Visible = 1;
                    app.DensityValueEditFieldLabel.Visible = 1;
                    app.RollingBallSizeSpinner.Visible = 0;
                    app.RollingBallSizeSpinnerLabel.Visible = 0;
            end
            box = str2num(app.BoxSelectionDropDown.Value);
            app.background_token(box) = 0;
            app.single_box_callback = 1;
            UpdateDisplay(app)
            app.single_box_callback = 0;
        end

        % Value changed function: DensityValueEditField
        function DensityValueEditFieldValueChanged(app, event)
            value = app.DensityValueEditField.Value;
            app.background_token = 0;
            app.single_box_callback = 1;
            UpdateDisplay(app)
            app.single_box_callback = 0;
        end

        % Value changed function: MedianFilterSizeSpinner
        function MedianFilterSizeSpinnerValueChanged(app, event)
            value = app.MedianFilterSizeSpinner.Value;
            box = str2num(app.BoxSelectionDropDown.Value);
            app.filtered_inset(box) = 0;
            app.single_box_callback = 1;
            UpdateDisplay(app)
            app.single_box_callback = 0;
        end

        % Value changed function: ApplyFilterCheckBox
        function ApplyFilterCheckBoxValueChanged(app, event)
            value = app.ApplyFilterCheckBox.Value;
            if value
                app.MedianFilterSizeSpinner.Enable = 'on';
                app.MedianFilterSizeSpinnerLabel.Enable = 'on';
            else
                app.MedianFilterSizeSpinner.Enable = 'off';
                app.MedianFilterSizeSpinnerLabel.Enable = 'off';
            end
            app.single_box_callback = 1;
            UpdateDisplay(app)
            app.single_box_callback = 0;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create GelBoxUIFigure and hide until all components are created
            app.GelBoxUIFigure = uifigure('Visible', 'off');
            app.GelBoxUIFigure.Colormap = [0.2431 0.149 0.6588;0.2431 0.1529 0.6745;0.2471 0.1569 0.6863;0.2471 0.1608 0.698;0.251 0.1647 0.7059;0.251 0.1686 0.7176;0.2549 0.1725 0.7294;0.2549 0.1765 0.7412;0.2588 0.1804 0.749;0.2588 0.1804 0.7608;0.2627 0.1882 0.7725;0.2588 0.1882 0.7804;0.2627 0.1961 0.7922;0.2667 0.2 0.8039;0.2667 0.2039 0.8157;0.2706 0.2078 0.8235;0.2706 0.2157 0.8353;0.2706 0.2196 0.8431;0.2745 0.2235 0.851;0.2745 0.2275 0.8627;0.2745 0.2314 0.8706;0.2745 0.2392 0.8784;0.2784 0.2431 0.8824;0.2784 0.2471 0.8902;0.2784 0.2549 0.898;0.2784 0.2588 0.902;0.2784 0.2667 0.9098;0.2784 0.2706 0.9137;0.2784 0.2745 0.9216;0.2824 0.2824 0.9255;0.2824 0.2863 0.9294;0.2824 0.2941 0.9333;0.2824 0.298 0.9412;0.2824 0.3059 0.9451;0.2824 0.3098 0.949;0.2824 0.3137 0.9529;0.2824 0.3216 0.9569;0.2824 0.3255 0.9608;0.2824 0.3294 0.9647;0.2784 0.3373 0.9686;0.2784 0.3412 0.9686;0.2784 0.349 0.9725;0.2784 0.3529 0.9765;0.2784 0.3569 0.9804;0.2784 0.3647 0.9804;0.2745 0.3686 0.9843;0.2745 0.3765 0.9843;0.2745 0.3804 0.9882;0.2706 0.3843 0.9882;0.2706 0.3922 0.9922;0.2667 0.3961 0.9922;0.2627 0.4039 0.9922;0.2627 0.4078 0.9961;0.2588 0.4157 0.9961;0.2549 0.4196 0.9961;0.251 0.4275 0.9961;0.2471 0.4314 1;0.2431 0.4392 1;0.2353 0.4431 1;0.2314 0.451 1;0.2235 0.4549 1;0.2196 0.4627 0.9961;0.2118 0.4667 0.9961;0.2078 0.4745 0.9922;0.2 0.4784 0.9922;0.1961 0.4863 0.9882;0.1922 0.4902 0.9882;0.1882 0.498 0.9843;0.1843 0.502 0.9804;0.1843 0.5098 0.9804;0.1804 0.5137 0.9765;0.1804 0.5176 0.9725;0.1804 0.5255 0.9725;0.1804 0.5294 0.9686;0.1765 0.5333 0.9647;0.1765 0.5412 0.9608;0.1765 0.5451 0.9569;0.1765 0.549 0.9529;0.1765 0.5569 0.949;0.1725 0.5608 0.9451;0.1725 0.5647 0.9412;0.1686 0.5686 0.9373;0.1647 0.5765 0.9333;0.1608 0.5804 0.9294;0.1569 0.5843 0.9255;0.1529 0.5922 0.9216;0.1529 0.5961 0.9176;0.149 0.6 0.9137;0.149 0.6039 0.9098;0.1451 0.6078 0.9098;0.1451 0.6118 0.9059;0.1412 0.6196 0.902;0.1412 0.6235 0.898;0.1373 0.6275 0.898;0.1373 0.6314 0.8941;0.1333 0.6353 0.8941;0.1294 0.6392 0.8902;0.1255 0.6471 0.8902;0.1216 0.651 0.8863;0.1176 0.6549 0.8824;0.1137 0.6588 0.8824;0.1137 0.6627 0.8784;0.1098 0.6667 0.8745;0.1059 0.6706 0.8706;0.102 0.6745 0.8667;0.098 0.6784 0.8627;0.0902 0.6824 0.8549;0.0863 0.6863 0.851;0.0784 0.6902 0.8471;0.0706 0.6941 0.8392;0.0627 0.698 0.8353;0.0549 0.702 0.8314;0.0431 0.702 0.8235;0.0314 0.7059 0.8196;0.0235 0.7098 0.8118;0.0157 0.7137 0.8078;0.0078 0.7176 0.8;0.0039 0.7176 0.7922;0 0.7216 0.7882;0 0.7255 0.7804;0 0.7294 0.7765;0.0039 0.7294 0.7686;0.0078 0.7333 0.7608;0.0157 0.7333 0.7569;0.0235 0.7373 0.749;0.0353 0.7412 0.7412;0.051 0.7412 0.7373;0.0627 0.7451 0.7294;0.0784 0.7451 0.7216;0.0902 0.749 0.7137;0.102 0.7529 0.7098;0.1137 0.7529 0.702;0.1255 0.7569 0.6941;0.1373 0.7569 0.6863;0.1451 0.7608 0.6824;0.1529 0.7608 0.6745;0.1608 0.7647 0.6667;0.1686 0.7647 0.6588;0.1725 0.7686 0.651;0.1804 0.7686 0.6471;0.1843 0.7725 0.6392;0.1922 0.7725 0.6314;0.1961 0.7765 0.6235;0.2 0.7804 0.6157;0.2078 0.7804 0.6078;0.2118 0.7843 0.6;0.2196 0.7843 0.5882;0.2235 0.7882 0.5804;0.2314 0.7882 0.5725;0.2392 0.7922 0.5647;0.251 0.7922 0.5529;0.2588 0.7922 0.5451;0.2706 0.7961 0.5373;0.2824 0.7961 0.5255;0.2941 0.7961 0.5176;0.3059 0.8 0.5059;0.3176 0.8 0.498;0.3294 0.8 0.4863;0.3412 0.8 0.4784;0.3529 0.8 0.4667;0.3686 0.8039 0.4549;0.3804 0.8039 0.4471;0.3922 0.8039 0.4353;0.4039 0.8039 0.4235;0.4196 0.8039 0.4118;0.4314 0.8039 0.4;0.4471 0.8039 0.3922;0.4627 0.8 0.3804;0.4745 0.8 0.3686;0.4902 0.8 0.3569;0.5059 0.8 0.349;0.5176 0.8 0.3373;0.5333 0.7961 0.3255;0.5451 0.7961 0.3176;0.5608 0.7961 0.3059;0.5765 0.7922 0.2941;0.5882 0.7922 0.2824;0.6039 0.7882 0.2745;0.6157 0.7882 0.2627;0.6314 0.7843 0.251;0.6431 0.7843 0.2431;0.6549 0.7804 0.2314;0.6706 0.7804 0.2235;0.6824 0.7765 0.2157;0.698 0.7765 0.2078;0.7098 0.7725 0.2;0.7216 0.7686 0.1922;0.7333 0.7686 0.1843;0.7451 0.7647 0.1765;0.7608 0.7647 0.1725;0.7725 0.7608 0.1647;0.7843 0.7569 0.1608;0.7961 0.7569 0.1569;0.8078 0.7529 0.1529;0.8157 0.749 0.1529;0.8275 0.749 0.1529;0.8392 0.7451 0.1529;0.851 0.7451 0.1569;0.8588 0.7412 0.1569;0.8706 0.7373 0.1608;0.8824 0.7373 0.1647;0.8902 0.7373 0.1686;0.902 0.7333 0.1765;0.9098 0.7333 0.1804;0.9176 0.7294 0.1882;0.9255 0.7294 0.1961;0.9373 0.7294 0.2078;0.9451 0.7294 0.2157;0.9529 0.7294 0.2235;0.9608 0.7294 0.2314;0.9686 0.7294 0.2392;0.9765 0.7294 0.2431;0.9843 0.7333 0.2431;0.9882 0.7373 0.2431;0.9961 0.7412 0.2392;0.9961 0.7451 0.2353;0.9961 0.7529 0.2314;0.9961 0.7569 0.2275;0.9961 0.7608 0.2235;0.9961 0.7686 0.2196;0.9961 0.7725 0.2157;0.9961 0.7804 0.2078;0.9961 0.7843 0.2039;0.9961 0.7922 0.2;0.9922 0.7961 0.1961;0.9922 0.8039 0.1922;0.9922 0.8078 0.1922;0.9882 0.8157 0.1882;0.9843 0.8235 0.1843;0.9843 0.8275 0.1804;0.9804 0.8353 0.1804;0.9765 0.8392 0.1765;0.9765 0.8471 0.1725;0.9725 0.851 0.1686;0.9686 0.8588 0.1647;0.9686 0.8667 0.1647;0.9647 0.8706 0.1608;0.9647 0.8784 0.1569;0.9608 0.8824 0.1569;0.9608 0.8902 0.1529;0.9608 0.898 0.149;0.9608 0.902 0.149;0.9608 0.9098 0.1451;0.9608 0.9137 0.1412;0.9608 0.9216 0.1373;0.9608 0.9255 0.1333;0.9608 0.9333 0.1294;0.9647 0.9373 0.1255;0.9647 0.9451 0.1216;0.9647 0.949 0.1176;0.9686 0.9569 0.1098;0.9686 0.9608 0.1059;0.9725 0.9686 0.102;0.9725 0.9725 0.0941;0.9765 0.9765 0.0863;0.9765 0.9843 0.0824];
            app.GelBoxUIFigure.Position = [100 100 1702 611];
            app.GelBoxUIFigure.Name = 'GelBox';
            app.GelBoxUIFigure.CloseRequestFcn = createCallbackFcn(app, @GelBoxUIFigureCloseRequest, true);

            % Create FileMenu
            app.FileMenu = uimenu(app.GelBoxUIFigure);
            app.FileMenu.Text = 'File';

            % Create LoadImageMenu
            app.LoadImageMenu = uimenu(app.FileMenu);
            app.LoadImageMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadImageButtonPushed, true);
            app.LoadImageMenu.Text = 'Load Image';

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
            app.OpticalDensitiesPanel.Position = [848 289 848 314];

            % Create raw_density
            app.raw_density = uiaxes(app.OpticalDensitiesPanel);
            xlabel(app.raw_density, 'Optical Density (A.U.)')
            ylabel(app.raw_density, 'Pixel')
            app.raw_density.Position = [163 47 254 209];

            % Create background_corrected_raw_density
            app.background_corrected_raw_density = uiaxes(app.OpticalDensitiesPanel);
            xlabel(app.background_corrected_raw_density, 'Optical Density (A.U.)')
            ylabel(app.background_corrected_raw_density, 'Pixel')
            app.background_corrected_raw_density.YColor = [0.9412 0.9412 0.9412];
            app.background_corrected_raw_density.Position = [401 47 254 209];

            % Create box_inset
            app.box_inset = uiaxes(app.OpticalDensitiesPanel);
            xlabel(app.box_inset, 'Optical Density')
            ylabel(app.box_inset, 'Pixel')
            app.box_inset.XColor = [0.9412 0.9412 0.9412];
            app.box_inset.YColor = [0.9412 0.9412 0.9412];
            app.box_inset.Position = [4 47 150 209];

            % Create BoxZoomLabel
            app.BoxZoomLabel = uilabel(app.OpticalDensitiesPanel);
            app.BoxZoomLabel.HorizontalAlignment = 'center';
            app.BoxZoomLabel.WordWrap = 'on';
            app.BoxZoomLabel.Position = [34 260 126 22];
            app.BoxZoomLabel.Text = 'Box Zoom';

            % Create RawOpticalDensityLabel_2
            app.RawOpticalDensityLabel_2 = uilabel(app.OpticalDensitiesPanel);
            app.RawOpticalDensityLabel_2.HorizontalAlignment = 'center';
            app.RawOpticalDensityLabel_2.WordWrap = 'on';
            app.RawOpticalDensityLabel_2.Position = [238 260 126 22];
            app.RawOpticalDensityLabel_2.Text = 'Raw Optical Density';

            % Create BackgroundCorrectedOpticalDensityLabel_2
            app.BackgroundCorrectedOpticalDensityLabel_2 = uilabel(app.OpticalDensitiesPanel);
            app.BackgroundCorrectedOpticalDensityLabel_2.HorizontalAlignment = 'center';
            app.BackgroundCorrectedOpticalDensityLabel_2.WordWrap = 'on';
            app.BackgroundCorrectedOpticalDensityLabel_2.Position = [477 257 126 28];
            app.BackgroundCorrectedOpticalDensityLabel_2.Text = 'Background Corrected Optical Density';

            % Create TotalAreaEditFieldLabel
            app.TotalAreaEditFieldLabel = uilabel(app.OpticalDensitiesPanel);
            app.TotalAreaEditFieldLabel.Position = [665 241 59 22];
            app.TotalAreaEditFieldLabel.Text = 'Total Area';

            % Create TotalAreaField
            app.TotalAreaField = uieditfield(app.OpticalDensitiesPanel, 'numeric');
            app.TotalAreaField.ValueDisplayFormat = '%.2f';
            app.TotalAreaField.Editable = 'off';
            app.TotalAreaField.HorizontalAlignment = 'center';
            app.TotalAreaField.Position = [743 241 100 22];

            % Create RollingBallSizeSpinnerLabel
            app.RollingBallSizeSpinnerLabel = uilabel(app.OpticalDensitiesPanel);
            app.RollingBallSizeSpinnerLabel.WordWrap = 'on';
            app.RollingBallSizeSpinnerLabel.Position = [665 153 52 33];
            app.RollingBallSizeSpinnerLabel.Text = 'Rolling Ball Size';

            % Create RollingBallSizeSpinner
            app.RollingBallSizeSpinner = uispinner(app.OpticalDensitiesPanel);
            app.RollingBallSizeSpinner.Limits = [1 Inf];
            app.RollingBallSizeSpinner.ValueChangedFcn = createCallbackFcn(app, @RollingBallSizeSpinnerValueChanged, true);
            app.RollingBallSizeSpinner.Position = [743 157 100 22];
            app.RollingBallSizeSpinner.Value = 80;

            % Create BackgroundCorrAreaLabel
            app.BackgroundCorrAreaLabel = uilabel(app.OpticalDensitiesPanel);
            app.BackgroundCorrAreaLabel.WordWrap = 'on';
            app.BackgroundCorrAreaLabel.Position = [665 65 81 30];
            app.BackgroundCorrAreaLabel.Text = 'Background Corr. Area';

            % Create BackgroundCorrAreaField
            app.BackgroundCorrAreaField = uieditfield(app.OpticalDensitiesPanel, 'numeric');
            app.BackgroundCorrAreaField.ValueDisplayFormat = '%.2f';
            app.BackgroundCorrAreaField.Editable = 'off';
            app.BackgroundCorrAreaField.HorizontalAlignment = 'center';
            app.BackgroundCorrAreaField.Position = [743 69 100 22];

            % Create BackgroundAreaLabel
            app.BackgroundAreaLabel = uilabel(app.OpticalDensitiesPanel);
            app.BackgroundAreaLabel.WordWrap = 'on';
            app.BackgroundAreaLabel.Position = [665 104 81 30];
            app.BackgroundAreaLabel.Text = 'Background Area';

            % Create BackgroundAreaField
            app.BackgroundAreaField = uieditfield(app.OpticalDensitiesPanel, 'numeric');
            app.BackgroundAreaField.ValueDisplayFormat = '%.2f';
            app.BackgroundAreaField.Editable = 'off';
            app.BackgroundAreaField.HorizontalAlignment = 'center';
            app.BackgroundAreaField.Position = [743 108 100 22];

            % Create DensityValueEditFieldLabel
            app.DensityValueEditFieldLabel = uilabel(app.OpticalDensitiesPanel);
            app.DensityValueEditFieldLabel.WordWrap = 'on';
            app.DensityValueEditFieldLabel.Visible = 'off';
            app.DensityValueEditFieldLabel.Position = [665 153 62 30];
            app.DensityValueEditFieldLabel.Text = 'Density Value';

            % Create DensityValueEditField
            app.DensityValueEditField = uieditfield(app.OpticalDensitiesPanel, 'numeric');
            app.DensityValueEditField.Limits = [0 Inf];
            app.DensityValueEditField.ValueDisplayFormat = '%.2f';
            app.DensityValueEditField.ValueChangedFcn = createCallbackFcn(app, @DensityValueEditFieldValueChanged, true);
            app.DensityValueEditField.HorizontalAlignment = 'center';
            app.DensityValueEditField.Visible = 'off';
            app.DensityValueEditField.Position = [743 156 100 22];
            app.DensityValueEditField.Value = 10;

            % Create BackgroundSubtractionDropDownLabel
            app.BackgroundSubtractionDropDownLabel = uilabel(app.OpticalDensitiesPanel);
            app.BackgroundSubtractionDropDownLabel.WordWrap = 'on';
            app.BackgroundSubtractionDropDownLabel.Position = [665 193 114 37];
            app.BackgroundSubtractionDropDownLabel.Text = 'Background Subtraction';

            % Create BackgroundSubtractionDropDown
            app.BackgroundSubtractionDropDown = uidropdown(app.OpticalDensitiesPanel);
            app.BackgroundSubtractionDropDown.Items = {'Rolling Ball', 'Linear', 'Constant Value', 'Cubic Spline', 'Polynomial', 'Polyline'};
            app.BackgroundSubtractionDropDown.ValueChangedFcn = createCallbackFcn(app, @BackgroundSubtractionDropDownValueChanged, true);
            app.BackgroundSubtractionDropDown.Position = [743 201 100 22];
            app.BackgroundSubtractionDropDown.Value = 'Rolling Ball';

            % Create ApplyFilterCheckBox
            app.ApplyFilterCheckBox = uicheckbox(app.OpticalDensitiesPanel);
            app.ApplyFilterCheckBox.ValueChangedFcn = createCallbackFcn(app, @ApplyFilterCheckBoxValueChanged, true);
            app.ApplyFilterCheckBox.Text = 'Apply Filter';
            app.ApplyFilterCheckBox.Position = [56 48 82 22];

            % Create MedianFilterSizeSpinnerLabel
            app.MedianFilterSizeSpinnerLabel = uilabel(app.OpticalDensitiesPanel);
            app.MedianFilterSizeSpinnerLabel.HorizontalAlignment = 'right';
            app.MedianFilterSizeSpinnerLabel.Position = [15 17 102 22];
            app.MedianFilterSizeSpinnerLabel.Text = 'Median Filter Size';

            % Create MedianFilterSizeSpinner
            app.MedianFilterSizeSpinner = uispinner(app.OpticalDensitiesPanel);
            app.MedianFilterSizeSpinner.Limits = [3 Inf];
            app.MedianFilterSizeSpinner.ValueChangedFcn = createCallbackFcn(app, @MedianFilterSizeSpinnerValueChanged, true);
            app.MedianFilterSizeSpinner.Position = [124 17 55 22];
            app.MedianFilterSizeSpinner.Value = 3;

            % Create GelImagePanel
            app.GelImagePanel = uipanel(app.GelBoxUIFigure);
            app.GelImagePanel.Title = 'Gel Image';
            app.GelImagePanel.Position = [6 6 833 597];

            % Create gel_image_axis
            app.gel_image_axis = uiaxes(app.GelImagePanel);
            app.gel_image_axis.XTick = [];
            app.gel_image_axis.YTick = [];
            app.gel_image_axis.Box = 'on';
            app.gel_image_axis.Position = [12 14 810 516];

            % Create AdjustImageButton
            app.AdjustImageButton = uibutton(app.GelImagePanel, 'push');
            app.AdjustImageButton.ButtonPushedFcn = createCallbackFcn(app, @AdjustImageButtonPushed, true);
            app.AdjustImageButton.Position = [16 537 100 22];
            app.AdjustImageButton.Text = 'Adjust Image';

            % Create NewBoxButton
            app.NewBoxButton = uibutton(app.GelImagePanel, 'push');
            app.NewBoxButton.ButtonPushedFcn = createCallbackFcn(app, @NewBoxButtonPushed, true);
            app.NewBoxButton.Position = [130 537 100 22];
            app.NewBoxButton.Text = 'New Box';

            % Create DeleteBoxButton
            app.DeleteBoxButton = uibutton(app.GelImagePanel, 'push');
            app.DeleteBoxButton.ButtonPushedFcn = createCallbackFcn(app, @DeleteBoxButtonPushed, true);
            app.DeleteBoxButton.Enable = 'off';
            app.DeleteBoxButton.Position = [247 537 101 22];
            app.DeleteBoxButton.Text = 'Delete Box';

            % Create BoxSelectionDropDownLabel
            app.BoxSelectionDropDownLabel = uilabel(app.GelImagePanel);
            app.BoxSelectionDropDownLabel.HorizontalAlignment = 'center';
            app.BoxSelectionDropDownLabel.Position = [354 537 98 22];
            app.BoxSelectionDropDownLabel.Text = 'Box Selection';

            % Create BoxSelectionDropDown
            app.BoxSelectionDropDown = uidropdown(app.GelImagePanel);
            app.BoxSelectionDropDown.Items = {};
            app.BoxSelectionDropDown.ValueChangedFcn = createCallbackFcn(app, @BoxSelectionDropDownValueChanged, true);
            app.BoxSelectionDropDown.Placeholder = 'No Data';
            app.BoxSelectionDropDown.Position = [451 537 100 22];
            app.BoxSelectionDropDown.Value = {};

            % Create FittingPanel
            app.FittingPanel = uipanel(app.GelBoxUIFigure);
            app.FittingPanel.Title = 'Fitting';
            app.FittingPanel.Position = [848 6 848 284];

            % Create background_corrected_raw_density_fit
            app.background_corrected_raw_density_fit = uiaxes(app.FittingPanel);
            xlabel(app.background_corrected_raw_density_fit, 'Optical Density (A.U.)')
            ylabel(app.background_corrected_raw_density_fit, 'Pixel')
            app.background_corrected_raw_density_fit.YColor = [0.9412 0.9412 0.9412];
            app.background_corrected_raw_density_fit.Position = [255 13 254 209];

            % Create raw_density_fit
            app.raw_density_fit = uiaxes(app.FittingPanel);
            xlabel(app.raw_density_fit, 'Optical Density (A.U.)')
            ylabel(app.raw_density_fit, 'Pixel')
            app.raw_density_fit.Position = [11 13 254 209];

            % Create RawOpticalDensityLabel
            app.RawOpticalDensityLabel = uilabel(app.FittingPanel);
            app.RawOpticalDensityLabel.HorizontalAlignment = 'center';
            app.RawOpticalDensityLabel.WordWrap = 'on';
            app.RawOpticalDensityLabel.Position = [88 221 126 22];
            app.RawOpticalDensityLabel.Text = 'Raw Optical Density';

            % Create BackgroundCorrectedOpticalDensityLabel
            app.BackgroundCorrectedOpticalDensityLabel = uilabel(app.FittingPanel);
            app.BackgroundCorrectedOpticalDensityLabel.HorizontalAlignment = 'center';
            app.BackgroundCorrectedOpticalDensityLabel.WordWrap = 'on';
            app.BackgroundCorrectedOpticalDensityLabel.Position = [338 218 126 28];
            app.BackgroundCorrectedOpticalDensityLabel.Text = 'Background Corrected Optical Density';

            % Create DrawFittingCheckBox
            app.DrawFittingCheckBox = uicheckbox(app.FittingPanel);
            app.DrawFittingCheckBox.ValueChangedFcn = createCallbackFcn(app, @DrawFittingCheckBoxValueChanged, true);
            app.DrawFittingCheckBox.Text = 'Draw Fitting';
            app.DrawFittingCheckBox.Position = [710 181 86 22];

            % Create FittingParametersButton
            app.FittingParametersButton = uibutton(app.FittingPanel, 'push');
            app.FittingParametersButton.ButtonPushedFcn = createCallbackFcn(app, @FittingParametersButtonPushed, true);
            app.FittingParametersButton.Position = [699 223 113 23];
            app.FittingParametersButton.Text = 'Fitting Parameters';

            % Create rsquaredField
            app.rsquaredField = uieditfield(app.FittingPanel, 'numeric');
            app.rsquaredField.ValueDisplayFormat = '%.3f';
            app.rsquaredField.Editable = 'off';
            app.rsquaredField.HorizontalAlignment = 'center';
            app.rsquaredField.Position = [626 181 70 22];

            % Create RsquaredLabel
            app.RsquaredLabel = uilabel(app.FittingPanel);
            app.RsquaredLabel.WordWrap = 'on';
            app.RsquaredLabel.Position = [550 181 72 22];
            app.RsquaredLabel.Text = 'R - squared';

            % Create NumberofBandsSpinnerLabel
            app.NumberofBandsSpinnerLabel = uilabel(app.FittingPanel);
            app.NumberofBandsSpinnerLabel.HorizontalAlignment = 'right';
            app.NumberofBandsSpinnerLabel.Position = [540 222 99 22];
            app.NumberofBandsSpinnerLabel.Text = 'Number of Bands';

            % Create NumberofBandsSpinner
            app.NumberofBandsSpinner = uispinner(app.FittingPanel);
            app.NumberofBandsSpinner.Limits = [1 Inf];
            app.NumberofBandsSpinner.ValueChangedFcn = createCallbackFcn(app, @NumberofBandsSpinnerValueChanged, true);
            app.NumberofBandsSpinner.Position = [647 222 46 24];
            app.NumberofBandsSpinner.Value = 1;

            % Create BandTable
            app.BandTable = uitable(app.FittingPanel);
            app.BandTable.ColumnName = {'Band No'; 'Color'; 'Area'; 'Relative Area'};
            app.BandTable.RowName = {};
            app.BandTable.Position = [518 20 321 148];

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