clearvars;
close all

r = 2;
d = 0.4;
T = 9/24;

L = 50;
dw = T/32;
L = ceil(L/dw)*dw;
w = 0:dw:L;
nw = length(w);
t_end = 10;
nt = 2^17;

h = w(2)-w(1);
t = linspace(0,t_end,nt);
dt = t(2)-t(1);
u0 = exp(-10*w);

u_int = w>=T;
w_int = w(u_int);

u = zeros(nw,nt);
u(:,1) = u0;

%%
A = zeros(nw-2);
A = A + diag(-d*ones(1,nw-2) - r*(w(2:end-1)>=T)) + diag(-1/(2*h)*ones(1,nw-3),1) + diag(1/(2*h)*ones(1,nw-3),-1);
A = sparse(A);
for i = 2:nt
    u(2:end-1,i) = u(2:end-1,i-1) + dt*A*u(2:end-1,i-1);
    u(2,i) = u(2,i) + dt*(1/(2*h))*u(1,i-1);
    u(end-1,i) = u(end-1,i) - dt*(1/(2*h))*u(end,i-1);
    u(1,i) = 2*r*trapz(w_int,u(u_int,i-1));
    u(end,i) = 0;
end

%%
figure;
ax = gca;
yL = min(u(:));
yR = max(u(:));
for i = 1:2^10:nt
    plot(ax,w,u(:,i)/trapz(w,u(:,i)))
%     ylim([yL,yR])
    drawnow
end