function [y_bands, y_fit,r_squared]= fit_gaussian(x,y,y_back,no_of_bands)

peaks=find_peaks('x',x, ...
    'y',y, ...
    'min_rel_delta_y',0.05, ...
    'min_x_index_spacing',2);

if numel(peaks.max_indices) == no_of_bands

    upper_band_x_estimate=peaks.max_indices(1);
    lower_band_x_estimate=peaks.max_indices(2);
else
    upper_band_x_estimate = 0.5*length(x);
    lower_band_x_estimate = 0.6*length(x);
end
target = y';

target = target - y_back;
[max_value,~]=max(target);


half_distance=(0.1*length(x));
alfa_estimate = -log(0.5)/(half_distance^2);

upper_curve_shape_estimate = alfa_estimate;
upper_amp_estimate=max_value;
lower_amp_estimate=upper_amp_estimate;
upper_skew_estimate=1;


par=[upper_band_x_estimate ...
    lower_band_x_estimate ...
    upper_curve_shape_estimate ...
    upper_amp_estimate  ...
    lower_amp_estimate ...
    upper_skew_estimate ...
    ];

j = 1;
e = [];
fig_disp = [];

opts=optimset('fminsearch');
opts.Display='off';
opts.MaxIter=1000;
opts.MaxFunEvals=10000;

[p_result,fval,exitflag,output] = fminsearch(@profile_error, par, opts);
p_result
r_squared = calculate_r_squared(y',y_fit+y_back);
if ~isempty(y_fit) && ~isempty(fig_disp) 
    figure(25)
    cla
    plot(y,x,'k','LineWidth',1.5)
    hold on
    plot(y(peaks.max_indices),peaks.max_indices,'bo')
    plot(y_fit+y_back,x,'gd')
    plot(y_back,x,'md')
    t = sprintf('r^2 = %.3f',r_squared);
    title(t)
    area_1 = trapz(x,y_bands(1,:));
    area_2 = trapz(x,y_bands(2,:));
    plot(y_bands(1,:)+y_back,x,'ro');
    plot(y_bands(2,:)+y_back,x,'bo');
    str1 = sprintf('Red: %.3f\n Blue: %.3f', area_1, area_2);

    xL=xlim;
    yL=ylim;
    text(0.99*xL(2),0.99*yL(2),str1,'HorizontalAlignment','right','VerticalAlignment','top')
end



    function trial_e = profile_error(par)

        [y_bands,y_fit] = calculate_profile(x,par);

        e(j) = 0;

        for i  = 1 : numel(target)
            e(j) = e(j) + (y_fit(i) - target(i))^2;
        end

        trial_e = e(end);

        if any(par<0)
            trial_e = 10^6;
        end

        j = j + 1;
        if fig_disp
            figure(22)
            cla
            plot(y_fit,x,'LineWidth',2)
            hold on
            plot(target,x,'LineWidth',2,'LineStyle','-.','Color','k')
            area_1 = trapz(x,y_bands(1,:));
            area_2 = trapz(x,y_bands(2,:));
            plot(y_bands(1,:),x,'ro');
            plot(y_bands(2,:),x,'bo');
            str1 = sprintf('Red: %.3f\n Blue: %.3f', area_1, area_2);

            xL=xlim;
            yL=ylim;
            text(0.99*xL(2),0.99*yL(2),str1,'HorizontalAlignment','right','VerticalAlignment','top')
        end
    end
    function [y_bands,y_fit] = calculate_profile(x,par)
        x1=par(1);
        x2=par(2);
        curve_shape1=par(3);
        amp1=par(4);
        amp2=par(5);
        skew1=par(6);
%         relative_shape=par(7);

        yred=skewed_Gaussian(x,x1,curve_shape1,amp1,skew1);
        yblue=skewed_Gaussian(x,x2,curve_shape1,amp2,skew1);

        y_fit=yblue+yred;
        y_bands(1,:) = yred;
        y_bands(2,:) = yblue;

    end
    function y=skewed_Gaussian(x,x0,gamma,A,skew1)
        offset=zeros(1,length(x));
        offset((x-x0)>0)=skew1*(x((x-x0)>0)-x0);

        y=A*exp(-gamma*(((x-x0)+offset).^2));
    end
end