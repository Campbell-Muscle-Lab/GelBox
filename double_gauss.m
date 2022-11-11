function [y1,y2,y_fit]= double_gauss(x,y,y_back)

fprintf('Double Gaussian optimization started....\n')

peaks=find_peaks('x',x, ...
    'y',y, ...
    'min_rel_delta_y',0.05, ...
    'min_x_index_spacing',5);

target = y';

target = target - y_back;
[max_value,max_index]=max(target);

half1= find(target(1:max_index)<(0.5*max_value),1,'last');
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

alfa_estimate = -log(0.5)/(half_distance^2);
skew_estimate = 1;
%findpeaks(data.y)

A1_guess = target(peaks.max_indices(1));
A2_guess = target(peaks.max_indices(2));
alfa_guess = alfa_estimate;
skew_guess = skew_estimate;
j = 1;
e = [];
par = [peaks.max_indices(1),peaks.max_indices(2),A1_guess,A2_guess,alfa_guess,skew_guess];


[p_result,fval,exitflag,output] = fminsearch(@profile_error, par, []);


    function trial_e = profile_error(par)


        [y1,y2,y_fit] = calculate_profile(x,par);

        e(j) = 0;

        for i  = 1 : numel(target)
            e(j) = e(j) + (y_fit(i) - target(i))^2;
        end

%         figure(3)
%         clf
% 
%         subplot(2,1,1)
%         plot(log10(e),'-r','MarkerFaceColor','r')
%         xlabel('Iteration')
%         ylabel('Log Error')

        trial_e = e(end);

        j = j + 1;

%         subplot(2,1,2)
%         plot(y_fit,'LineWidth',2)
%         hold on
%         plot(target,'LineWidth',2,'LineStyle','-.','Color','k')
%         plot(y1)
%         plot(y2)
    end

    function [y1,y2,y_fit] = calculate_profile(x,par)

        x1 = par(1);
        x2 = par(2);
        A1 = par(3);
        A2 = par(4);
        alfa = par(5);
        skew = par(6);

        y1 = gaussian_fun(x,x1,alfa,A1,skew);
        y2 = gaussian_fun(x,x2,1.5*alfa,A2,skew);

        y_fit = y1 + y2;

    end

    function y = gaussian_fun(x,x0,alfa,A,skew)

        offset=zeros(1,length(x));
        offset((x-x0)>0)=skew*(x((x-x0)>0)-x0);
        y=A*exp(-alfa*(((x-x0)+offset).^2));

    end
end