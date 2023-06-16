% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;

addpath("../../../ODEFittingFns/")

make_save = true;

opts.force_serial = false;
opts.n_starts = 1;
opts.raw_error_opts.assume_independent_time_series = true;

cohort_name = "cohort_230124175743017";


[p,lb,ub] = basePars();
npars = numel(p);

fn = @computeTimeSeries;
fn_opts.condition_on_previous = false;

weights = 1;
optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%% load ABM data
files.data = sprintf("../../data/%s/summary_short.mat",cohort_name);

P = optimizeSMParsFromABM(files,p,fn,fn_opts,lb,ub,optim_opts,weights,opts);

if make_save
    save("data/SMFitToABM_New.mat","P","cohort_name")
end

%% reset path
rmpath("../../../ODEFittingFns/")
