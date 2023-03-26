function dx = odefn(x,p,death_rate)

% the ODE function defining the SM

dx = [-1;1] * p(1) * x(1) + p(2) * [2 - sum(x)/p(3);-1] * x(2) - death_rate.*x;
