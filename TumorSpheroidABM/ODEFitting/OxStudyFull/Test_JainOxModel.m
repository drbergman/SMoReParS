clearvars;

lambda = 5.96;
alphaRP = 0.99;
theta = 1;
VT = 5.48;
V0 = .55;
alphaP = 8.69;
kalpha = 21.67;
a = 0.5;
rho0 = 0.06;
delta0 = 0.96;
kdelta = 9.22;
b = 1;
alphaR = 3.0;
p = [lambda;alphaRP;theta;VT;V0;alphaP;kalpha;a;rho0;delta0;kdelta;b;alphaR];

IC = [.10;.90;0;0];
dose = 0;

tt = linspace(0,3,100);
out = computeTimeSeries_jain(p,tt,dose);

figure;
plot(tt,out(:,1))
yyaxis right
plot(tt,out(:,2))
ylim([0 1])
legend({"Total","Proliferating Proportion"},"location","northwest")
% set(gca,'YScale','log')