function cm_pars = drawCMParameters(x,par_names,D,T)

cm_pars = cell(length(x),1);
for i = 1:length(x)
    if ~any(par_names(i) == D.keys)
        continue;
    end
    cm_pars{i} = icdf(D(par_names(i)),x(i));
    if any(par_names(i) == T.keys)
        cm_pars{i} = feval(T(par_names(i)),cm_pars{i});
    end
end

end