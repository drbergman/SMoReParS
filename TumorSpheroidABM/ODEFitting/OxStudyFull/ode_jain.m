function dx = ode_jain(x,p,dose)

% dose = oxaliplatin concentration
% p = [lambda,alphaRP,theta,VT,V0,alphaP,kalpha,a,rho0,delta0,kdelta,b,alphaR]
VF = p(4)-sum(x);
dx = zeros(4,1);

if VF<=0
    disp('')
end
%     reentry = 0;
% else
    reentry = p(2)*x(2)/(1+(p(5)/VF)^p(3));
% end


dx(1:2) = dx(1:2) + p(1)*x(1)*[-1;2] + ... % proliferation
    reentry*[1;-1];  % reentry
dx([1,3]) = dx([1,3]) + p(6)*x(1)/(1+(p(7)/dose)^p(8))*[-1;1] + ... % arrest of proliferating cells
    p(9)*x(3)*[1;-1]; % reactivation into proliferating compartment (this should perhaps instead go into resting compartment)
dx([2,4]) = dx([2,4]) + p(13)*x(2)/(1+(p(7)/dose)^p(8))*[-1;1] + ... % arrest of resting cells
    p(9)*x(4)*[1;-1]; % reactivation into resting compartment
dx(3:4) = dx(3:4) - p(10)*x(3:4)/(1+(p(11)/dose)^p(12)); % apoptosis of arrested cells


