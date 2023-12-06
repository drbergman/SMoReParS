function out = sampleFromSM(x,par_names,BS,T,D,vals,nsamps,fn,C,fn_opts,sum_fn)

% runs the SM on a LHS of ODE parameter space as defined by x.
% par_names stores the correspondence between
% indices in x and parameters in the SM. D is a vector of distributions for the
% parameters, so that a call to ICDF can return the appropriate value. T is
% a dictionary of transformation maps for the parameters. e.g. converting
% to an integer.

n_abm_pars = length(par_names);
n_sm_pars = size(BS,1);
abm_pars = cell(n_abm_pars,1);
for i = 1:n_abm_pars
    transform_map = T(par_names(i));
    abm_pars{i} = transform_map(icdf(D(par_names(i)),x(i)));
end

colons = repmat({':'},n_abm_pars,1);
Vq = zeros(n_sm_pars,2);
for pi = 1:n_sm_pars
    Vq(pi,1) = interpn(vals{:},squeeze(BS(pi,colons{:},1)),abm_pars{:});
    Vq(pi,2) = interpn(vals{:},squeeze(BS(pi,colons{:},2)),abm_pars{:});
end

X = lhsdesign(nsamps,n_sm_pars,"Smooth","off"); % LHS on [0,1]
points = (1-X).*Vq(:,1)' + X.*Vq(:,2)'; % use these as weights to interpolate between the boundaries

out = zeros(nsamps,1);

for i = 1:nsamps
    out(i) = fn(points(i,:)',C,fn_opts);
end

out = sum_fn(out);

