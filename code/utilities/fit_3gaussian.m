function [y_bands, y_fit,r_squared]= fit_3gaussian(gel_data,x,y,y_back,box_no)
                        
            par_est = struct();
            par_con = struct();
            
            box_pars = struct();
            box_pars = gel_data.par_est(box_no);
            
            box_cons = struct();
            box_cons = gel_data.par_con(box_no);
            
            names = fieldnames(box_pars);
            for m = 1 : numel(names)
                par_est.(names{m}) = box_pars.(names{m});
                par_con.(names{m}) = box_cons.(names{m});

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
            
            shape_and_skew = 0;
            shape = 0;
            skew = 0;
            
            if (numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1) || ...
                    (any(par_con.shape_parameter) && ...
                    any(par_con.skew_parameter))
                par(9) = par_est.shape_parameter(2);
                par(10) = par_est.skew_parameter(2);
                par(11) = par_est.shape_parameter(3);
                par(12) = par_est.skew_parameter(3);
                shape_and_skew = 1;
            elseif (numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) == 1) || ...
                    (any(par_con.shape_parameter) && ...
                    ~any(par_con.skew_parameter))
                par(9) = par_est.shape_parameter(2);
                par(10) = par_est.shape_parameter(3);
                shape = 1;
            elseif (numel(unique(par_est.shape_parameter)) == 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1) || ...
                    ~(any(par_con.shape_parameter) && ...
                    any(par_con.skew_parameter))
                par(9) = par_est.skew_parameter(2);
                par(10) = par_est.skew_parameter(3);
                skew = 1;
            end

            
            j = 1;
            e = [];
            

            target = y';

            target = target - y_back;

            no_of_parameters = numel(par);
            
            A_constraints = zeros(1,no_of_parameters);
            A_constraints(1) = 1;
            B_constants = [0];
            
            lower_bounds=zeros(1,no_of_parameters);
            upper_bounds=Inf*ones(1,no_of_parameters);
            
            no_of_bands = numel(par_con.(names{m}));
            
            l = 1;
            for u = 1 : 1
                for m = 2 : numel(names)
                    if par_con.(names{m})(u)
                        lower_bounds(l) = par_est.(names{m})(u);
                        upper_bounds(l) = par_est.(names{m})(u);
                    end
                    l = l + 1;
                end
            end
            
            if par_con.peak_location(2)
                lower_bounds(5) = par_est.peak_location(2);
                upper_bounds(5) = par_est.peak_location(2);
            end
            
            if par_con.amplitude(2)
                lower_bounds(6) = par_est.amplitude(2);
                upper_bounds(6) = par_est.amplitude(2);
            end
            
            if par_con.peak_location(3)
                lower_bounds(7) = par_est.peak_location(3);
                upper_bounds(7) = par_est.peak_location(3);
            end
            
            if par_con.amplitude(3)
                lower_bounds(8) = par_est.amplitude(3);
                upper_bounds(8) = par_est.amplitude(3);
            end
            
            if shape_and_skew
                if par_con.shape_parameter(2)
                    lower_bounds(9) = par_est.shape_parameter(2);
                    upper_bounds(9) = par_est.shape_parameter(2);
                end
                if par_con.skew_parameter(2)
                    lower_bounds(10) = par_est.skew_parameter(2);
                    upper_bounds(10) = par_est.skew_parameter(2);
                end
                if par_con.shape_parameter(3)
                    lower_bounds(11) = par_est.shape_parameter(2);
                    upper_bounds(11) = par_est.shape_parameter(2);
                end
                if par_con.skew_parameter(3)
                    lower_bounds(12) = par_est.skew_parameter(2);
                    upper_bounds(11) = par_est.skew_parameter(2);
                end
            end
            
            if shape
                if par_con.shape_parameter(2)
                    lower_bounds(9) = par_est.shape_parameter(2);
                    upper_bounds(9) = par_est.shape_parameter(2);
                end
                if par_con.shape_parameter(3)
                    lower_bounds(10) = par_est.shape_parameter(3);
                    upper_bounds(10) = par_est.shape_parameter(3);
                end
            end

            if skew
                if par_con.skew_parameter(2)
                    lower_bounds(9) = par_est.skew_parameter(2);
                    upper_bounds(9) = par_est.skew_parameter(2);
                end

                if par_con.skew_parameter(3)
                    lower_bounds(10) = par_est.skew_parameter(2);
                    upper_bounds(10) = par_est.skew_parameter(2);
                end
            end
            

            opts=optimset('fminsearch');
            opts.Display='off';
            opts.MaxIter=1000;
            opts.MaxFunEvals=10000;
            
            [p_result,fval,exitflag,output]= ...
                fminsearchcon(@profile_error_3gaussian,par, ...
                lower_bounds,upper_bounds,[],[], ...
                [],opts);
            
            par_fit.band_no = [1;2;3];
            par_fit.peak_location(1,1) = p_result(1);
            par_fit.shape_parameter(1,1) = p_result(2);
            par_fit.amplitude(1,1) = p_result(3);
            par_fit.skew_parameter(1,1) = p_result(4);
            par_fit.peak_location(2,1) = p_result(5);
            par_fit.amplitude(2,1) = p_result(6);
            par_fit.peak_location(3,1) = p_result(7);
            par_fit.amplitude(3,1) = p_result(8);
            
            if (numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1) || ...
                    (any(par_con.shape_parameter) && ...
                    any(par_con.skew_parameter))
                par_fit.shape_parameter(2,1) = p_result(9);
                par_fit.skew_parameter(2,1) = p_result(10);
                par_fit.shape_parameter(3,1) = p_result(11);
                par_fit.skew_parameter(3,1) = p_result(12);
            elseif (numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) == 1) || ...
                    (any(par_con.shape_parameter) && ...
                    ~any(par_con.skew_parameter))
                par_fit.shape_parameter(2,1) = p_result(9);
                par_fit.skew_parameter(2,1) = par_fit.skew_parameter(1,1);
                par_fit.shape_parameter(3,1) = p_result(10);
                par_fit.skew_parameter(3,1) = par_fit.skew_parameter(1,1);
            elseif (numel(unique(par_est.shape_parameter)) == 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1) || ...
                    ~(any(par_con.shape_parameter) && ...
                    any(par_con.skew_parameter))
                par_fit.shape_parameter(2,1) = par_fit.shape_parameter(1,1);
                par_fit.skew_parameter(2,1) = p_result(9);
                par_fit.shape_parameter(3,1) = par_fit.shape_parameter(1,1);
                par_fit.skew_parameter(3,1) = p_result(10);
            else
                par_fit.shape_parameter(2,1) = par_fit.shape_parameter(1,1);
                par_fit.skew_parameter(2,1) = par_fit.skew_parameter(1,1);
                par_fit.shape_parameter(3,1) = par_fit.shape_parameter(1,1);
                par_fit.skew_parameter(3,1) = par_fit.skew_parameter(1,1);
                
            end

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

                if (numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1) || ...
                    (any(par_con.shape_parameter) && ...
                    any(par_con.skew_parameter))
                    curve_shape2 = par(9);
                    skew2 = par(10);
                    curve_shape3 = par(11);
                    skew3 = par(12);
                elseif (numel(unique(par_est.shape_parameter)) ~= 1 && ...
                    numel(unique(par_est.skew_parameter)) == 1) || ...
                    (any(par_con.shape_parameter) && ...
                    ~any(par_con.skew_parameter))
                    curve_shape2 = par(9);
                    skew2 = skew1;
                    curve_shape3 = par(10);
                    skew3 = skew1;
                elseif (numel(unique(par_est.shape_parameter)) == 1 && ...
                    numel(unique(par_est.skew_parameter)) ~= 1) || ...
                    ~(any(par_con.shape_parameter) && ...
                    any(par_con.skew_parameter))
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