function out = computeTimeSeries_jain(p,tt,dose,~)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [52;38;8;2] so that [90;10] in the two G1/S and G2/M groups. this split was chosen by weighting the first by [11,8] and the second by [4,1], the durations of the four phases;

sol = ode45(@(t,x) ode_jain(x,p,dose),[0 3],[.1;.9;0;0]);
temp = deval(sol,tt)';
total = sum(temp,2);
out = [total,sum(temp(:,[1,3]),2)./total];
