function out = computeTimeSeries(p,tt,chemo_conc,chemo_death_is_continuous)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [90;10];
if chemo_death_is_continuous
    death_rate = p(4) * chemo_conc;
else
    death_rate = p(4) * chemo_conc * [p(1);p(2)];
end
sol = ode45(@(t,x) odefn(x,p,death_rate),[0 3],[90;10]);

out = deval(sol,tt)';