function [y_bands, y_fit]= fit_gaussian(x,y,y_back,no_of_bands)

peaks=find_peaks('x',x, ...
    'y',y, ...
    'min_rel_delta_y',0.01, ...
    'min_x_index_spacing',2);

if numel(peaks.max_indices) ~= no_of_bands
    warn_text = sprintf('The number of peaks is not equaled to number of bands.\nThe number of bands is adjusted to the number of peaks');
    warndlg(warn_text)
    no_of_bands = numel(peaks.max_indices);
end
target = y';

target = target - y_back;
[max_value,max_index]=max(target);

par = zeros(no_of_bands,3);

half_distance=(0.05*length(x));
alfa_estimate = -log(0.5)/(half_distance^2);
skew = 1;

for m = 1:no_of_bands
par(m,1) = target(peaks.max_indices(m));
par(m,2) = alfa_estimate;
par(m,3) = peaks.max_indices(m);
par(m,4) = skew;
end

par = reshape(par.',1,[]);

j = 1;
e = [];

[p_result,fval,exitflag,output] = fminsearch(@profile_error, par, []);

    function trial_e = profile_error(par)

        [y_bands,y_fit] = calculate_profile(x,par);

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
        
        if any(par<0)
            trial_e = 10^6;
        end

        j = j + 1;

%         subplot(2,1,2)
%         plot(y_fit,'LineWidth',2)
%         hold on
%         plot(target,'LineWidth',2,'LineStyle','-.','Color','k')
%         for i = 1:no_of_bands
%         plot(y_bands(i,:));
%         end
    end
    function [y_bands,y_fit] = calculate_profile(x,par)       
        y_fit = 0;
        k = 1;
        for i = 1 : no_of_bands
            par_dum(i,:) = par(k:k+3);
            k = k + numel(par_dum(i,:));
       end

        for i = 1 : no_of_bands

            offset=zeros(1,length(x));
            offset((x-par_dum(i,3))>0)=par_dum(i,4)*(x((x-par_dum(i,3))>0)-par_dum(i,3));
            y_bands(i,:) = par_dum(i,1)*exp(-par_dum(i,2)*(((x-par_dum(i,3))+offset).^2));
            y_fit = y_fit + y_bands(i,:);
        end

    end
end