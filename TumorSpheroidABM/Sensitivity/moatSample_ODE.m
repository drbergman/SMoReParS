function out = moatSample_ODE(x,par_names,D)

% runs the ODE at the point in parameter space defined by x.
% par_names stores the correspondence between
% indices in x and parameters in M. D is a vector of distributions for the
% parameters, so that a call to ICDF can return the appropriate value.
for i = 1:numel(x)
    v_temp = icdf(D(par_names(i)),x(i));
    switch par_names(i)
        case "lambda"
            p(1) = v_temp;
        case "alpha"
            p(2) = v_temp;
        case "K"
            p(3) = v_temp;
        case "delta"
            p(4) = v_temp;
        case "g1_prop0"
            p(5) = v_temp;
        otherwise
            error("Have not yet planned for %s to be varied.",par_names(i))
    end
end

M.save_pars.make_save = false; % make sure I'm not creating save files for any of this
M.plot_pars.plot_fig = false; % make sure I'm not plotting anything either

out = sum(computeTimeSeries(p,3));


