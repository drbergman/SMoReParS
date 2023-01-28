clearvars;

% p.lambda = 1;
% p.alpha = 1;
% p.K = 1e3;
% p.delta = .1;
% p.g1_prop0 = 0.8;
% p.P = @(x) 1 - sum(x)/p.K;
p(1) = 1;
p(2) = 1;
p(3) = 1e3;
p(4) = .1;
p(5) = 0.4;

tt = linspace(0,3,5000);
out = computeTimeSeries(p,tt);

sol = ode45(@(t,x) odefn(x,p),[0 3],100*[p(5);1-p(5)]);

figure;
subplot(2,1,1)
plot(sol.x,[sol.y;sum(sol.y)])
subplot(2,1,2)
plot(tt,out)