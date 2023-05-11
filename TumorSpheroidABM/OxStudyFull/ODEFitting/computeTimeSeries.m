function out = computeTimeSeries(p,tt,dose,opts,~)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [90;10];

if opts.link_phase_death_rates
    if ~opts.phase_dependent_death
        death_rate = p(4) * dose;
    else
        death_rate = p(4) * dose * [p(1);p(2)];
    end
else
    if opts.hill_activation
        death_rate = p(4:5).*dose^p(6)/(dose^p(6)+p(7)^p(6));
    else
        death_rate = p(4:5) * dose;
    end
end

sol2 = ode45(@(t,x) odefn(x,p,death_rate),[0 tt(end)],[90;10]);
sol = ode23(@(t,x) odefn(x,p,death_rate),[0 tt(end)],[90;10]);
temp = deval(sol,tt)';
total = sum(temp,2);
out = [total,temp(:,2)./total];
