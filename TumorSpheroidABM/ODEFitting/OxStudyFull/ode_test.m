function dx = ode_test(x,p,T,fn)

% p = [lambda,alphaRP,theta,VT,V0,alphaP,kalpha,a,rho0,delta0,kdelta,b]
VF = p(4)-sum(x);
dx = zeros(3,1);
reactivation = p(2)*x(2)*fn(VF,p(5),p(3),p(4));
arrest = p(6)/((p(7)/T)^p(8))*x(1);
dx(1) = -p(1)*x(1) + reactivation - arrest + p(9)*x(3);
dx(2) = 2*p(1)*x(1) - reactivation;
dx(3) = arrest - p(10)/((p(11)/T)^p(12)+1) - p(9)*x(3);
