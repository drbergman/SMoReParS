% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;

addpath("../../../ODEFittingFns/")

make_save = true;

force_serial = false;
n_starts = 1;
raw_error_opts = struct();
raw_error_opts.assume_independent_time_series = true;

cohort_name = "cohort_2401151523";


[p,lb,ub] = basePars();
npars = numel(p);

sm.fn = @computeTimeSeries;
sm.opts.condition_on_previous = false;

weights = 1;
optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%% load ABM data
files.data = sprintf("../../data/%s/summary.mat",cohort_name);

P = optimizeSMParsFromABM(files,sm,p,lb,ub,optim_opts,weights,...
    force_serial=force_serial,n_starts=n_starts,raw_error_opts=raw_error_opts);

if make_save
    save("data/SMFitToABM_New.mat","P","cohort_name")
end

%% reset path
rmpath("../../../ODEFittingFns/")
