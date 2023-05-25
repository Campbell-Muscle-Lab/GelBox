function [y_bands, y_fit,r_squared]= fit_gaussian2(x,y,y_back,no_of_bands)

target  = [];
target = y';

target = target - y_back;

[max_value,max_index] = max(target);

half1= ...
    find(target(1:max_index)<(0.5*max_value),1,'last');
% Second one
half2=max_index-1 + ...
    find(target(max_index:end)<(0.5*max_value), ...
    1,'first');
% Half-distance
if ((length(half1)>0)&&(length(half2)>0))
    half_distance= ...
        min([abs(max_index-half1) abs(max_index-half2)]);
else
    % Guess something plausible
    half_distance=(0.1*length(image_data.zoomed_pixels));
end
%buraya diger half distance

alfa_estimate = -log(0.5)/(half_distance^2);

amp_estimate = max_value;

skew_estimate = 1;

rel_shape_estimate = 1.5;

peaks=find_peaks('x',x, ...
    'y',target, ...
    'min_rel_delta_y',0.01, ...
    'min_x_index_spacing',1);

if numel(peaks.max_indices) ~= no_of_bands
    if no_of_bands == 2
        peaks.max_indices = [];
        peaks.max_indices(1) = round(0.6*length(x));
        peaks.max_indices(2) = round(0.5*length(x));
    end
end

band_estimate_1 = peaks.max_indices(1);
band_estimate_2 = peaks.max_indices(2);
shape_estimate = alfa_estimate;
amp_1 = max_value;
amp_2 = amp_1;
skew = skew_estimate;
relative_shape_estimate = rel_shape_estimate;

 
par =[band_estimate_1...
    band_estimate_2...
    shape_estimate...
    amp_1...
    amp_2...
    skew...
    relative_shape_estimate];

A_constraints=[1 -1 0 0 0 0 0];
B_constants=[0];
lower_bounds=zeros(1,length(par));
upper_bounds=Inf*ones(1,length(par));

upper_bounds(1)=length(x);
upper_bounds(2)=length(x);
upper_bounds(4)=max(target);
 upper_bounds(5)=max(target);

lower_bounds(7)=0.05;
upper_bounds(7)=5;
j = 1;
e = [];
opts=optimset('fminsearch')
opts.Display='off';
opts.MaxIter=1000;
opts.MaxFunEvals=10000;
[p_result,fval,exitflag,output]= ...
                fminsearchcon(@profile_error,par, ...
                lower_bounds,upper_bounds,[],[],[],opts);
fprintf('utku')
r_squared = calculate_r_squared(y,y_fit+y_back);
if ~isempty(y_fit)
    figure(2)
    cla
    plot(y,x,'k','LineWidth',1.5)
    hold on
    plot(y(peaks.max_indices),peaks.max_indices,'bo')
    plot(y_fit+y_back,x,'rd')
    plot(y_back,x,'md')
    t = sprintf('r^2 = %.3f',r_squared);
    title(t)
end
    function trial_e = profile_error(par)

        [y_1,y_2,y_fit] = calculate_profile(x,par);

        e(j) = 0;

        for i  = 1 : numel(target)
            e(j) = e(j) + (y_fit(i) - target(i))^2;
        end

%         figure(33)
%         clf
% 
%         subplot(2,1,1)
%         plot(log10(e),'-r','MarkerFaceColor','r')
%         xlabel('Iteration')
%         ylabel('Log Error')

        trial_e = e(end);
        j = j+1

%         subplot(2,1,2)
%         plot(y_fit,'LineWidth',2)
%         hold on
%         plot(target,'LineWidth',2,'LineStyle','-.','Color','k')
%         plot(y_1);
%         plot(y_2)
    y_bands(1,:) = y_1;
    y_bands(2,:) = y_2;

    end

    function [y_1,y_2,y_fit] = calculate_profile(x,par)

        x_1 = par(1);
        x_2 = par(2);
        curve_shape_1 = par(3);
        amplitude_1 = par(4);
        amplitude_2 = par(5);
        skew_1 = par(6);
        relative_shape = par(7);

        y_1 = skewed_gauss(x, x_1, curve_shape_1, amplitude_1,skew_1);
        y_2 = skewed_gauss(x, x_2, relative_shape*curve_shape_1, amplitude_2,skew_1);

        y_fit = y_1 + y_2;
    end

    function s_g = skewed_gauss(x,x0,alfa,A,skew1)

        offset=zeros(1,length(x));
        offset((x-x0)>0)=skew1*(x((x-x0)>0)-x0);

        s_g=A*exp(-alfa*(((x-x0)+offset).^2));

    end

% par = zeros(no_of_bands,3);
%
% for m = 1:no_of_bands
%     par(m,1) = max_value;
%     par(m,2) = alfa_estimate;
%     par(m,3) = peaks.max_indices(m);
%     par(m,4) = skew;
% end
%
% par = reshape(par.',1,[]);
%
% j = 1;
% e = [];
%
% [p_result,fval,exitflag,output] = fminsearch(@profile_error, par, []);
%
% r_squared = calculate_r_squared(y,y_fit+y_back);
% if ~isempty(y_fit)
%     figure(2)
%     cla
%     plot(y,x,'k','LineWidth',1.5)
%     hold on
%     plot(y(peaks.max_indices),peaks.max_indices,'bo')
%     %     plot(y_fit+y_back,x,'rd')
%     plot(y_back,x,'md')
%     t = sprintf('r^2 = %.3f',r_squared);
%     title(t)
% end
%
%
%
%     function trial_e = profile_error(par)
%
%         [y_bands,y_fit] = calculate_profile(x,par);
%
%         e(j) = 0;
%
%         for i  = 1 : numel(target)
%             e(j) = e(j) + (y_fit(i) - target(i))^2;
%         end
%
%         trial_e = e(end);
%
%         if any(par<0)
%             trial_e = 10^6;
%         end
%
%         j = j + 1;
%     end

end