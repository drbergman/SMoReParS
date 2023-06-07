function out = computeTimeSeries(p,tt,dose,fn_opts,~)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [90;10;0;0];

if isfield(fn_opts,"p_setup_fn")
    p = fn_opts.p_setup_fn(p);
end

if dose==0
    dose_arrest_factor = 0;
    dose_apoptosis_factor = 0;
    recovery_rate = 0;
else
    dose_arrest_factor = 1/(1+(p(6)/dose)^p(7));
    if numel(p)==11
        dose_apoptosis_factor = p(8)/(1+(p(9)/dose)^p(10));
        recovery_rate = p(11);
    else
        switch dose
            case 0.75
                dose_apoptosis_factor = p(8);
            case 7.55
                dose_apoptosis_factor = p(8) + p(9);
            otherwise
                error("Did not expect this dose")
        end
        recovery_rate = p(10);
    end
end

sol = ode45(@(t,x) odefn(x,p,dose_arrest_factor,dose_apoptosis_factor,recovery_rate),[0 tt(end)],[90;10;0;0]);
temp = deval(sol,tt)';
total = sum(temp,2);
out = [total,(temp(:,2)+temp(:,4))./total];
