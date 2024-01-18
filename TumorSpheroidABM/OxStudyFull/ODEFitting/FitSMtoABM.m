% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;

file_name = "SMFitToABM_LMS_bounded";
make_save = false;
addpath("../../../SurrogateModelFns/")
addpath("~/Documents/MATLAB/myfunctions/")

force_serial = true;
n_starts = 1;
checkpoint_filename = "data/temp_optimal";
save_every_iter = 20; % wait at least this many iterations between saves
save_every_sec = 2*60; % wait at least this many seconds between saves
% raw_error_opts.assume_independent_time_series = true; % assume that the two time series have diagonal covariance matrices at each time point

cohort_name = "cohort_2401151702";
files.data = sprintf("../../data/%s/summary.mat",cohort_name);
% files.previous_optim_file = "data/temp_optimal.mat";

load("data/SMFitToData_LMS_bounded.mat","fixed_pars","lb","ub","model_type","optim_opts")
load("data/SMFitToData_LMS_bounded.mat","fn","fn_opts","sm")
if ~exist("sm","var")
    sm.custom_solve_sm_fn = @(sm,p,t,c,d,~,~) computeTimeSeries(p,t,c,sm.opts,[]);
    sm.opts = fn_opts;
    % sm.fn = fn;
    % sm.opts = fn_opts;
end
optim_opts.Display = "off";
[p,~,~,~] = fixParameters(model_type,fixed_pars);

weights = [1;1;1];

%% optimize sm pars

P = optimizeSMParsFromABM(files,sm,p,lb,ub,optim_opts,weights,...
    force_serial=force_serial,n_starts=n_starts,checkpoint_filename=checkpoint_filename,...
    save_every_iter=save_every_iter,save_every_sec=save_every_sec);

if make_save
    save("data/" + file_name,"P","cohort_name") 
end

rmpath("../../../SurrogateModelFns/")
