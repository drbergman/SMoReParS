% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;

make_save = true;
addpath("../../../ODEFittingFns/")
addpath("~/Documents/MATLAB/myfunctions/")

model_type = "LogisticModel";


opts.force_serial = false;
opts.n_starts = 10;
opts.temp_profile_name = "data/temp_optimal";
opts.save_every_iter = 100; % wait at least this many iterations between saves
opts.save_every_sec = 10*60; % wait at least this many seconds between saves
opts.raw_error_opts.assume_independent_time_series = true; % assume that the two time series have diagonal covariance matrices at each time point

cohort_name = "cohort_2305311216";
files.data_file = sprintf("../../data/%s/summary.mat",cohort_name);
files.previous_optim_file = "data/temp_optimal.mat";

load("data/ODEFitToData.mat","fixed_pars")

fn = @computeTimeSeries;
[p,lb,ub,fn_opts.p_setup_fn] = fixParameters(model_type,fixed_pars);

optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

weights = [1;1;1];

%% optimize sm pars
P = optimizeSMParsFromABM(files,p,fn,fn_opts,lb,ub,optim_opts,weights,opts);

if make_save
    save("data/ODEFitToABM.mat","P") 
end

rmpath("../../../ODEFittingFns/")
