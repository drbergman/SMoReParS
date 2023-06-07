function [p,lb,ub,f,fixed_vals] = fixParameters(model_type,fixed_pars)

D = parameterOrdering(model_type);

[p,lb,ub] = basePars(model_type,fixed_pars);
npars = numel(p);
fixed_inds = zeros(numel(fixed_pars),1);
for i = 1:numel(fixed_pars)
    fixed_inds(i) = D(fixed_pars(i));
end
lb(fixed_inds) = [];
ub(fixed_inds) = [];
fixed_vals = p(fixed_inds);
p(fixed_inds) = [];
f = @(p) p_setup_fn(p,fixed_inds,fixed_vals,npars);

end

function p = p_setup_fn(p_in,fixed_inds,fixed_vals,npars)

p = NaN(npars,1);
p(fixed_inds) = fixed_vals;
p(isnan(p)) = p_in;

end
