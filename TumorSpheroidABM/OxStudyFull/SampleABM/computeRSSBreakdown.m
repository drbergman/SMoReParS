function rss = computeRSSBreakdown(E,D)

% E: experimental data as struct array; each element is for an experimental condition
%   fields:
%   A: mean time series
%   S: standard deviation time series
% D: simulated data as array of size [nt,n_conditions,n_time_series]

n_conditions = numel(E.C);
n_time_series = size(E.D(1).A,2);
n_time_points = size(E.D(1).A,1);
rss = nan(n_time_points,n_time_series,n_conditions);

for ci = 1:n_conditions
    % rss(:,ci) = sum(((E.D(ci).A-squeeze(D(:,ci,:)))./E.D(ci).S).^2,1);
    rss(:,:,ci) = ((E.D(ci).A-squeeze(D(:,ci,:)))./E.D(ci).S).^2;
end
