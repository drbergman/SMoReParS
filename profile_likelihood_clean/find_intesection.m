function [values] = find_intesection(array_x, array_y, value)
% This function finds the x-value of the intersection between a function 
% and a value. That is it finds, via interpolation, the value of x so that 
% f(x) = value.

N = length(array_x);
% number of intersection values
num_of_values = 0;

for i = 1:N-1
    if (array_y(i) - value)*(array_y(i+1) - value) < 0
        % linearly interpolate to find the intersection value of array_x
        slope = (array_y(i+1) - array_y(i)) / (array_x(i+1) - array_x(i));
        values(num_of_values + 1) = array_x(i) + (value - array_y(i)) / slope;

        % update number of intersection values
        num_of_values = num_of_values + 1;
    elseif (array_y(i) == value)
        values(num_of_values + 1) = array_x(i);

        % update number of intersection values
        num_of_values = num_of_values + 1;
    end
end

% if there was only one value that crossed the threshold, add last
% parameter value as the upper bound for the confidence interval
if (num_of_values == 1)
    values(2) = array_x(N);
end
% if no parameter values crossed the threshold, return an empty array
if (~exist("values", "var"))
    values = [];
end

end