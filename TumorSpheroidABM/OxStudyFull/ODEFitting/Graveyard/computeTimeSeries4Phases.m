function out = computeTimeSeries4Phases(p,tt,dose,opts)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [52;38;8;2] so that [90;10] in the two G1/S and G2/M groups. this split was chosen by weighting the first by [11,8] and the second by [4,1], the durations of the four phases;

if dose==0 % then effective rate is the transition rate
    eff_rate = [p(1);p(3)];
else

    if opts.link_phase_death_rates
        eff_rate = (1-dose*p(6))*[p(1);p(3)];
    else
        if opts.hill_activation
            eff_rate = (1-p(6:7)./(1+(p(9)/dose)^p(8))).*[p(1);p(3)];
        else
            eff_rate = (1-p(6:7).*dose).*[p(1);p(3)];
        end
    end
end
sol = ode45(@(t,x) odefn_4phases(x,p,eff_rate),[0 3],[52;38;8;2]);

temp = deval(sol,tt)';
total = sum(temp,2);
out = [total,sum(temp(:,3:4),2)./total];
