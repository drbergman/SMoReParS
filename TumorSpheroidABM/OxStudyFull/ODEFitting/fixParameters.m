function [p,lb,ub,f] = fixParameters(model_type,fixed_pars)

D = parameterOrdering(model_type);

[p,lb,ub] = basePars(fixed_pars);
fixed_inds = zeros(numel(fixed_pars),1);
for i = 1:numel(fixed_pars)
    fixed_inds(i) = D(fixed_pars(i));
end
lb(fixed_inds) = [];
ub(fixed_inds) = [];
fixed_vals = p(fixed_inds);
p(fixed_inds) = [];
f = @(p) p_setup_fn(p,fixed_inds,fixed_vals);

end

function p = p_setup_fn(p_in,fixed_inds,fixed_vals)

p = zeros(11,1);
p(fixed_inds) = fixed_vals;
p(p==0) = p_in;

end
