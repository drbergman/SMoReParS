function dx = odefn_4phases(x,p,eff_rate)

% the ODE function defining the SM with all four phases

% dx = [-1;1] * p(1) * x(1) + p(2) * [2 - sum(x)/p(3);-1] * x(2) - death_rate.*x;
dx = zeros(4,1);
dx(1) = -p(1)*x(1) + (2-sum(x)/p(5)) * p(4) * x(4);
dx(2) = eff_rate(1)*x(1) - p(2) * x(2); % eff_rate(1) = (1-C*p(6))*p(1)
dx(3) = p(2)*x(2) - p(3) * x(3);
dx(4) = eff_rate(2)*x(3) - p(4)*x(4); % eff_rate(1) = (1-C*p(6))*p(3)
