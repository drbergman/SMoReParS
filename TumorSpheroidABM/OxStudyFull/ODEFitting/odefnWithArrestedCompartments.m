function dx = odefnWithArrestedCompartments(x,p,dose_arrest_factor,dose_apoptosis_factor)

% p = [lambda,alphaRP,theta,VT,V0,alphaP,kalpha,a,rho0,delta0,kdelta,b,alphaR]
VF = p(4)-sum(x);
dx = zeros(4,1); % [G2/M; G1/S; arrested G2/M; arrested G1/S];

if VF<=0
    reentry = 0;
else
    reentry = p(2)*x(2)/(1+(p(5)/VF)^p(3));
end

dx(1:2) = dx(1:2) + p(1)*x(1)*[-1;2] + ... % proliferation
    reentry*[1;-1];  % reentry
dx([1,3]) = dx([1,3]) + p(6)*x(1)*dose_arrest_factor*[-1;1] + ... % arrest of proliferating cells
    p(9)*x(3)*[1;-1]; % reactivation into proliferating compartment (this should perhaps instead go into resting compartment)
dx([2,4]) = dx([2,4]) + p(13)*x(2)*dose_arrest_factor*[-1;1] + ... % arrest of resting cells
    p(9)*x(4)*[1;-1]; % reactivation into resting compartment
dx(3:4) = dx(3:4) - dose_apoptosis_factor*x(3:4); % apoptosis of arrested cells


