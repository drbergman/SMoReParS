function [I1,I2] = getProfileInds(vals,threshold)

% takes in profile vals and returns the bounds on the indices that keeps 
% the profile within the threshold
% profile: a 1xn array of output values
% threshold: chi2inv value that determines difference between max permitted
% value and the min of the second row in profile

[I1,I2,~] = findMaxCrosses(vals,threshold);

if isempty(I1)
    I1 = 1;
elseif numel(I1)==1 % the next index is the first to drop below the max value
    I1 = I1 + 1;
else
    error("Too many drops found.")
end

if isempty(I2)
    I2 = length(vals);
elseif numel(I2)==1 % then this is actually the index to use
else
    error("Too many rises found.")
end
