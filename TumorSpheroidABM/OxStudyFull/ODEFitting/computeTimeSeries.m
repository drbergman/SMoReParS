function out = computeTimeSeries(p,tt,dose,fn_opts,~)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [90;10;0;0];

if isfield(fn_opts,"p_setup_fn")
    p = fn_opts.p_setup_fn(p);
end

dose_arrest_factor = 1/(1+(p(6)/dose)^p(7));
dose_apoptosis_factor = p(8)/(1+(p(9)/dose)^p(10));

sol = ode45(@(t,x) odefn(x,p,dose_arrest_factor,dose_apoptosis_factor),[0 tt(end)],[90;10;0;0]);
temp = deval(sol,tt)';
total = sum(temp,2);
out = [total,(temp(:,2)+temp(:,4))./total];
