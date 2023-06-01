function out = computeTimeSeriesWithArrestedCompartments(p,tt,dose,~,~)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [10;90;0;0];

% dose = oxaliplatin concentration

% p = [lambda,alphaRP,theta,VT,V0,alphaP,kalpha,a,rho0,delta0,kdelta,b,alphaR]

dose_arrest_factor = 1/(1+(p(7)/dose)^p(8));
dose_apoptosis_factor = p(10)/(1+(p(11)/dose)^p(12));
sol = ode45(@(t,x) odefnWithArrestedCompartments(x,p,dose_arrest_factor,dose_apoptosis_factor),[0 tt(end)],[10;90;0;0]);
temp = deval(sol,tt)';
total = sum(temp,2);
out = [total,sum(temp(:,[1,3]),2)./total];
