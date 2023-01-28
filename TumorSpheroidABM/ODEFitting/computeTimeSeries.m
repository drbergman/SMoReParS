function out = computeTimeSeries(p,tt)

% sol = ode45(@(t,x) odefn(x,p),[0 3],100 * ([0;1] + [1;-1] * p(5)));
sol = ode45(@(t,x) odefn(x,p),[0 3],[90;10]);
out = deval(sol,tt)';