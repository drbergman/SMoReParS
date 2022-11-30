function y = vbmodel_fix(varied_p, x0, index,fixed_p, time)
%VBMODEL_FIX Summary of this function goes here
%   Detailed explanation goes here

% create array of parameters
n = length(varied_p) + 1;
p = zeros(n,1);

for i = 1:n
    if i < index
        p(i) = varied_p(i);
    elseif i == index
        p(i) = fixed_p;
    else
        p(i) = varied_p(i-1);
    end
end
options = odeset('AbsTol',1e-6,'RelTol',1e-8);
[~,y] = ode15s(@(t,x) vbmodel(x,p),time,x0,options); % solve ODE systems
end

