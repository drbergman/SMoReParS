function out = moatSample_SM(x,par_names,D,BS,pars,nsamps)

% runs the ODE on a LHS of ODE parameter space as defined by x.
% par_names stores the correspondence between
% indices in x and parameters in M. D is a vector of distributions for the
% parameters, so that a call to ICDF can return the appropriate value.
PG = cell(7,1);
for i = 1:numel(x)
    v_temp = icdf(D(par_names(i)),x(i));
    switch par_names(i)
        case "carrying_capacity"
            PG{1} = v_temp;
        case "g1_to_s"
            PG{4} = v_temp;
        case "s_to_g2"
            PG{5} = v_temp;
        case "g2_to_m"
            PG{6} = v_temp;
        case "m_to_g1"
            PG{7} = v_temp;
        case "arrest_prob_g1"
            M.cycle_pars.arrest_prob_g1 = v_temp;
        case "arrest_prob_g2"
            M.cycle_pars.arrest_prob_g2 = v_temp;
        case "apop_rate"
            M.pars.apop_rate = v_temp;
        case "move_rate_microns"
            PG{3} = v_temp;
        case "occmax_2d"
            PG{2} = min(7,floor(v_temp)); % make sure that this value cannot exceed 7 (otherwise cells will proliferate even when no space exists)

        otherwise
            error("Have not yet planned for %s to be varied.",par_names(i))
    end
end

Vq = zeros(3,2);
for pi = 1:3
    Vq(pi,1) = interpn(pars{:},squeeze(BS(pi,:,:,:,:,:,:,:,1)),PG{:});
    Vq(pi,2) = interpn(pars{:},squeeze(BS(pi,:,:,:,:,:,:,:,2)),PG{:});
end

points = zeros(nsamps,3);
for pi=1:3
    temp = linspace(Vq(pi,1),Vq(pi,2),nsamps+2);
    temp([1,end]) = []; % remove end points for LHS
    points(:,pi) = temp(randperm(length(temp)));
end

out = 0;

for i = 1:nsamps
    out = out + sum(computeTimeSeries(points(i,:)',3));
end

out = out / nsamps;

