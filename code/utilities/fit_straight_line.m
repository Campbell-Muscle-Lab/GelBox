function [grad,intercept,r,p,y_fit]=fit_straight_line(x_data,y_data);
% Functions fits straight line to data
% Uses the stastistics toolbox

% Some error checking
[x_rows,x_cols]=size(x_data);
if (x_rows<x_cols)
    x_data=x_data';
end
[y_rows,y_cols]=size(y_data);
if (y_rows<y_cols)
    y_data=y_data';
end
if ~isequal(size(x_data),size(y_data))
    error('X and Y data sets are different sizes (fit_linear function)');
end

% Create the x-matrix (allows for intercepts)
no_of_points=length(x_data);
x_matrix=[ones(no_of_points,1) x_data];

% Fit
[coefficients,intervals,residuals,residual_intervals,stats]=regress(y_data,x_matrix);
grad=coefficients(2);
intercept=coefficients(1);
r_matrix=corrcoef(x_data,y_data);
r=r_matrix(1,2);
p=stats(3);

% Calculate y_fit
y_fit=intercept+grad*x_data;
