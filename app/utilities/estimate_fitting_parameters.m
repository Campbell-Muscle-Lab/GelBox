function [par_est,par_con] = estimate_fitting_parameters(x,y,y_back, ...
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
                     
                     par_con.band_no  = 1;
                     par_con.peak_location = false;
                     par_con.shape_parameter = false;
                     par_con.amplitude = false;
                     par_con.skew_parameter = false;

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
                     
                     par_con.band_no  = [1;2];
                     par_con.peak_location = [false;false];
                     par_con.shape_parameter = [false;false];
                     par_con.amplitude = [false;false];
                     par_con.skew_parameter = [false;false];

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
                    
                     par_con.band_no  = [1;2;3];
                     par_con.peak_location = [false;false;false];
                     par_con.shape_parameter = [false;false;false];
                     par_con.amplitude = [false;false;false;];
                     par_con.skew_parameter = [false;false;false];
            end


            
        end