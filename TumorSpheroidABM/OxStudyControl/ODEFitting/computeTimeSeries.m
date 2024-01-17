function out = computeTimeSeries(p,tt,data, condition_on_previous, resample_t)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [90;10];

arguments
    p (:,1) double
    tt (:,1) double
    data double = [90;10] % as of writing this (24/01/15), I expect this value to be overwritten every time. the check for size(data,2)==1 below could cause issues if calls to this do not supply data
    condition_on_previous logical = false
    resample_t double = []
end

if condition_on_previous
    if size(data,1) < length(tt) - 1
        % then not enough data points are given
        error("Not enough data points given to condition on the previous data at each time point supplied.")
    end
    % out = zeros(length(tt),2);
    % out(1,:) = data(1,:);
    sol = ode45(@(t,x) odefn(x,p),tt(1:2),data(1,:)');
    % out(2,:) = deval(sol,tt(2));
    for i = 3:length(tt)
        sol = odextend(sol,[],tt(i),data(i-1,:)');
        % out(i,:) = deval(sol,tt(i));
    end
    if isempty(resample_t)
        out = deval(sol,tt*(1-eps()))';
    else
        out = deval(sol,resample_t*(1-eps()))';
    end
else
    if length(tt) == 1
        tspan = [0 tt];
    else
        tspan = tt;
    end
    sol = ode45(@(t,x) odefn(tt,x,p,{[]}),tspan,[90;10]); % do not put data here for [90;10] or else unexpected behavior could occur if data is supplied like I expect it is for some calls
    if isempty(resample_t)
        out = deval(sol,tt)';
    else
        out = deval(sol,resample_t)';
    end
end
if size(data,2)==1
    out = sum(out,2);
end

% sol = ode45(@(t,x) odefn(x,p),[0 3],[90;10]);
%
% out = deval(sol,tt)';