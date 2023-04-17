function [par_min,par_max] = getProfileBounds(profile,threshold)

% takes in profile and returns the bounds on the parameter that keeps the
% profile within the threshold
% profile: a 2xn array with row 1 being parameter values and row 2 being
% output values
% threshold: chi2inv value that determines difference between max permitted
% value and the min of the second row in profile

assert(size(profile,1)==2)

[I1,I2,max_val] = findMaxCrosses(profile(2,:),threshold);

if isempty(I1)
    par_min = profile(1,1);
elseif numel(I1)==1 % use linear interpolation to estimate where in this change the threshold is hit
    x = profile(1,I1+(0:1));
    y = profile(2,I1+(0:1));
    par_min = x(1) + diff(x)*(max_val-y(1))/diff(y);
else
    error("Too many drops found.")
end

if isempty(I2)
    par_max = profile(1,end);
elseif numel(I2)==1 % use linear interpolation to estimate where in this change the threshold is hit
    x = profile(1,I2+(0:1));
    y = profile(2,I2+(0:1));
    par_max = x(1) + diff(x)*(max_val-y(1))/diff(y);
else
    error("Too many rises found.")
end
