% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;

file_name = "SMFitToABM_LMS_bounded";
make_save = true;
addpath("../../../ODEFittingFns/")
addpath("~/Documents/MATLAB/myfunctions/")

opts.force_serial = false;
opts.n_starts = 10;
opts.temp_profile_name = "data/temp_optimal";
opts.save_every_iter = 20; % wait at least this many iterations between saves
opts.save_every_sec = 2*60; % wait at least this many seconds between saves
% opts.raw_error_opts.assume_independent_time_series = true; % assume that the two time series have diagonal covariance matrices at each time point

cohort_name = "cohort_2306062212";
files.data = sprintf("../../data/%s/summary.mat",cohort_name);
% files.previous_optim_file = "data/temp_optimal.mat";

load("data/SMFitToData_LMS_bounded.mat","fixed_pars","fn","lb","ub","fn_opts","model_type","optim_opts")
optim_opts.Display = "off";
[p,~,~,~] = fixParameters(model_type,fixed_pars);

weights = [1;1;1];

%% optimize sm pars
P = optimizeSMParsFromABM(files,p,fn,fn_opts,lb,ub,optim_opts,weights,opts);

if make_save
    save("data/" + file_name,"P","cohort_name") 
end

rmpath("../../../ODEFittingFns/")
