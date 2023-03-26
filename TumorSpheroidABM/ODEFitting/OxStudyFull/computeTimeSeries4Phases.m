function out = computeTimeSeries4Phases(p,tt,chemo_conc)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [52;38;8;2] so that [90;10] in the two G1/S and G2/M groups. this split was chosen by weighting the first by [11,8] and the second by [4,1], the durations of the four phases;
eff_rate = (1-chemo_conc*p(6))*[p(1);p(3)];
sol = ode45(@(t,x) odefn_4phases(x,p,eff_rate),[0 3],[52;38;8;2]);

out = deval(sol,tt)';