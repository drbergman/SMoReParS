function [inds_before_drop,inds_before_rise] = getProfileBoundaryInds(vals,threshold)

% takes in profile vals and returns the bounds on the indices that keeps 
% the profile within the threshold
% profile: a 1xn array of goodness-of-fit values
% threshold: chi2inv value that determines difference between max permitted
% value and the min of the second row in profile

[inds_before_drop,inds_before_rise,~] = findMaxCrosses(vals,threshold);

if isempty(inds_before_drop) % then it started below the threshold and never temporarily increased beyond the threshold
    inds_before_drop = 1;
elseif numel(inds_before_drop)==1 % the next index is the first to drop below the max value
    inds_before_drop = inds_before_drop + 1;
else
    error("Too many drops found.")
end

if isempty(inds_before_rise)
    inds_before_rise = length(vals);
elseif numel(inds_before_rise)==1 % then this is actually the index to use
else
    error("Too many rises found.")
end

if inds_before_drop > inds_before_rise
    error("The first drop happens after the first rise, indicating we start below the max")
end
