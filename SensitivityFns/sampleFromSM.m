% function out = sampleFromSM(x,par_names,BS,T,D,vals,nsamps,sm_functional,C,fn_opts,sum_fn)
function out = sampleFromSM(x, BS, vals, sm_functional, options)
% runs the SM on a LHS of ODE parameter space as defined by x.
% par_names stores the correspondence between
% indices in x and parameters in the SM. D is a vector of distributions for the
% parameters, so that a call to ICDF can return the appropriate value. T is
% a dictionary of transformation maps for the parameters. e.g. converting
% to an integer.

arguments
    x
    BS
    vals
    sm_functional
    options.par_names {mustBeText} = strings(1,numel(x))
    % options.T = dictionary()
    options.T = configureDictionary("string","function_handle")
    options.D = configureDictionary("string","prob.ProbabilityDistribution")
    options.nsamps {mustBeInteger} = 100
    options.sum_fn = @mean
end

n_abm_pars = length(options.par_names);
n_sm_pars = size(BS,1);
abm_pars = cell(n_abm_pars,1);
for i = 1:n_abm_pars
    if ~any(options.par_names(i) == options.D.keys)
        continue;
    end
    if any(options.par_names(i) == options.T.keys)
        transform_map = options.T(options.par_names(i));
    else
        transform_map = @(x) x;
    end
    abm_pars{i} = transform_map(icdf(options.D(options.par_names(i)),x(i)));
end

colons = repmat({':'},n_abm_pars,1);
Vq = zeros(n_sm_pars,2);
for pi = 1:n_sm_pars
    Vq(pi,1) = interpn(vals{:},squeeze(BS(pi,colons{:},1)),abm_pars{:});
    Vq(pi,2) = interpn(vals{:},squeeze(BS(pi,colons{:},2)),abm_pars{:});
end

X = lhsdesign(options.nsamps,n_sm_pars,"Smooth","off"); % LHS on [0,1]
points = (1-X).*Vq(:,1)' + X.*Vq(:,2)'; % use these as weights to interpolate between the boundaries

out = zeros(options.nsamps,1);

for i = 1:options.nsamps
    out(i) = sm_functional(points(i,:)');
end

out = options.sum_fn(out);