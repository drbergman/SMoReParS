function out = computeTimeSeries(p,tt)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [90;10];

% sol = ode45(@(t,x) odefn(x,p),[0 3],100 * ([0;1] + [1;-1] * p(5)));
sol = ode45(@(t,x) odefn(x,p),[0 3],[90;10]);
out = deval(sol,tt)';