% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;
opts.force_serial = false;
opts.raw_error_opts.assume_independent_time_series = false;

cohort_name = "cohort_230124175743017";

addpath("../../../ODEFittingFns/")

p = basePars();
npars = numel(p);

fn = @computeTimeSeries;
fn_opts.condition_on_previous = true;

lb = [0;0;0];
ub = [Inf;30;1e4];

weights = 1;
optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%% load ABM data
data_file = sprintf("../../data/%s/summary_short.mat",cohort_name);

P = optimizeSMParsFromABM(data_file,p,fn,fn_opts,lb,ub,optim_opts,weights,opts);

save("data/OptimalParameters_Using_DependentStatesAndTimeSeries.mat","P")

rmpath("../../../ODEFittingFns/")