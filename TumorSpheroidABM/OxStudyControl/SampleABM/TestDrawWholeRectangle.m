% A quick script to test drawing from an entire rectangle using a
% hyperplane to interpolate the log-likelihood

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

b1 = 2;
n = 101;
u = linspace(0,1,n);
n = length(u);
du = u(2)-u(1);

%% n = 1
A1 = -log(1+(exp(-b1)-1).*u)./b1;
figureOnRight;
plot(A1,u,"Marker",".")
axis([0 1 0 1])

%% n = 2
g2 = -20;
A2 = log(1+(exp(g2)-1).*u)./g2;
figureOnRight;
plot(A2,u,"Marker",".")
axis([0 1 0 1])
