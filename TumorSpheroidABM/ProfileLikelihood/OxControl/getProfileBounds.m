function [par_min,par_max] = getProfileBounds(profile,threshold)

% takes in profile and returns the bounds on the parameter that keeps the
% profile within the threshold
% profile: a 2xn array with row 1 being parameter values and row 2 being
% output values
% threshold: chi2inv value that determines difference between max permitted
% value and the min of the second row in profile

min_val = min(profile(2,:));
max_val = min_val+threshold;
changes = diff(profile(2,:)>max_val);

% I1 = find(changes==-1,1); % first time the profile dipped below the threshold, i.e. where the min should par value should be
I1 = find(changes==-1); % first time the profile dipped below the threshold, i.e. where the min should par value should be
if isempty(I1)
    par_min = profile(1,1);
elseif numel(I1)==1 % use linear interpolation to estimate where in this change the threshold is hit
    x = profile(1,I1+(0:1));
    y = profile(2,I1+(0:1));
    par_min = x(1) + diff(x)*(max_val-y(1))/diff(y);
else
    error("Too many drops found.")
end

% I2 = find(changes==1,1,"last"); % first time the profile rose above the threshold, i.e. where the max should par value should be
I2 = find(changes==1); % first time the profile rose above the threshold, i.e. where the max should par value should be
if isempty(I2)
    par_max = profile(1,end);
elseif numel(I2)==1 % use linear interpolation to estimate where in this change the threshold is hit
    x = profile(1,I2+(0:1));
    y = profile(2,I2+(0:1));
    par_max = x(1) + diff(x)*(max_val-y(1))/diff(y);
else
    error("Too many rises found.")
end
