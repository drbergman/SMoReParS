function LP = linkArrest(LP,cohort_pars)

if cohort_pars.link_arrest_coeffs % enforce that the arrest coefficients are identical
    for li = 1:numel(cohort_pars.linkings) % linking index
        arrest_coeff_inds = [];
        vals = {};
        paths = {};
        for i = 1:numel(LP)
            if ~iscell(LP(i).path) && startsWith(LP(i).path(end),"arrest_coeff") && any(endsWith(LP(i).path(end),cohort_pars.linkings{li}))
                arrest_coeff_inds(end+1) = i;
                vals{end+1} = LP(i).values;
                paths{end+1} = LP(i).path;
            end
        end
        if length(arrest_coeff_inds)>1 % only make this change if two arrest coeffs are varied
            assert(numel(unique(cellfun(@numel,vals)))==1) % make sure that each of these values has the same number of elements
            LP(end+1).values = cat(2,vals{:});
            LP(end).path = paths;
            LP(arrest_coeff_inds) = [];
        end
    end
end