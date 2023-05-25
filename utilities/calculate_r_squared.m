function r_sq = calculate_r_squared(y,y_fit)
% Given 2 arrays, y and y_fit, calculate the r_squared

[r1,c1]=size(y);
[r2,c2]=size(y_fit);

if (r1~=r2)
    y_fit=y_fit';
end

sum_y=sum(y);
sum_squares_residuals=sum((y_fit-y).^2);

y_mean=mean(y);

sum_squares_mean=sum((y-y_mean).^2);

r_sq = 1-sum_squares_residuals/sum_squares_mean;
