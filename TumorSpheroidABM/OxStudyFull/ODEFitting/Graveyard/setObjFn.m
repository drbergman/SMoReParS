function F = setObjFn(tt,data,data_std,C,objfn_opts)

if objfn_opts.is_abm_data % then data is from the ABM and we can compare with the SM directly
    F = @(p) sum(arrayfun(@(i) sum(((computeTimeSeries(p,tt,C(i),objfn_opts.phase_dependent_death)-data(:,:,i))./data_std(:,:,i)).^2,'all'),1:numel(C)),'all');
end