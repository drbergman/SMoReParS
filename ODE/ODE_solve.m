% this script solves the ODE, but it may be duplicated in
% TumorSpheroidABM/ODEFitting/TestODE.m

clearvars;

alpha21 = 0.16;
lambda = 0.5;
delta = 0.1;
K = 1e3;

p = @(x) 1-sum(x)/K;

F = @(x) [-lambda*x(1)+alpha21*(1+p(x))*x(2) - delta*x(1);...
           lambda*x(1)-alpha21*x(2)-delta*x(2)];

sol = ode45(@(t,x) F(x),[0 1000],[100;100]);

figure; plot(sol.x,sol.y)