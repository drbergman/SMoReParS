clearvars;

lambda = 9.57;
alphaRP = 20;
VT = 1.05e3;
V0 = 1.46e3;
theta = 9;
alphaP = 8.32;
kalpha = 8.10;
rho0 = 0.06;
delta0 = 1.33;
kdelta = 119.14;
a = 3;
b = 4;

IC = [0.9*1091*[.9;.8];0];
T = .005;
p = [lambda,alphaRP,theta,VT,V0,alphaP,kalpha,a,rho0,delta0,kdelta,b];

fn = @(VF,V0,theta,VT) 1/((V0/VF)^theta+1);
% fn = @(VF,V0,theta,VT) VF/VT;

sol = ode15s(@(t,x) ode_test(x,p,T,fn),[0 15],IC);

figure;
plot(sol.x,sol.y)
legend({"Proliferating","Resting","Arrested"})
% set(gca,'YScale','log')