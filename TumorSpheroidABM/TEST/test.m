clearvars;

w = 9/24;
r = 2;
f = @(x,y) [-1/w,2*r;1/w,-r]*[x;y];

sol = ode45(@(t,x) f(x(1),x(2)),[0 10],[1e8;1e5]);

figure;
plot(sol.x,sol.y(1,:)./sum(sol.y,1))

%%

a = -1/w;
b = r-1/w;
c = 2*r;

p = [a,b,c];
rh = roots(p);
rh = rh(rh>=0);

p = rh/(1+rh);
yline(p)