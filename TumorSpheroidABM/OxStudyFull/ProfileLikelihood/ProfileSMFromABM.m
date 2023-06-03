
% This script sets up and calls the profile likelihood method for the ODE
% SM on a grid of ABM parameters. compare_every allows for comparing with a
% subset of ABM time points. This could be improved slightly by
% interpolating the ABM output at those times rather than the rounding I do
% to get tt. The function profileLikelihood computes the profile likelihood
% of all ODE model parameters at a given ABM parameter vector.

clearvars;
addpath("../../../ODEFittingFns/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

addpath("~/Documents/MATLAB/myfunctions/")

file_name = "Profiles_SMFromABM_Fit_b";

% folder to store plots and text files
cohort_name = "cohort_2305311216";

files.par_file = "../ODEFitting/data/SMFitToABM_Fit_b.mat";
files.data_file = sprintf("../../data/%s/summary.mat",cohort_name);
files.previous_profile_file = "data/temp_profile.mat";
% files.previous_profile_file = "ProfileLikelihoods.mat";

options.profile_likelihood_options.save_all_pars = true;
options.force_serial = false;
options.temp_profile_name = "data/temp_profile";
options.save_every_iter = 50; % wait at least this many iterations between saves
options.save_every_sec = 5*60; % wait at least this many seconds between saves

load("../ODEFitting/data/SMFitToData_Fit_b.mat","fixed_pars","fn","lb","ub","fn_opts","model_type","optim_opts")
optim_opts.Display = "off";
[p,~,~,~] = fixParameters(model_type,fixed_pars);

%% setup profile params
n_sm_pars = numel(p);

profile_params.initial_step_prop = .01*ones(n_sm_pars,1);
profile_params.min_num_steps = 10*ones(n_sm_pars,1);
profile_params.smallest_par_step = 1e-1*ones(n_sm_pars,1); % do not let the step size go below this as it steps towards the boundary/threshold
profile_params.smallest_par_step(3) = 10; % use a larger min step for K
profile_params.shrinking_factor = 0.9; % factor by which to shrink dx as it gets close to enforced boundary
profile_params.threshold = chi2inv(0.95,n_sm_pars); % compute threshold value for the parameter confidence intervals

profile_params.secondary_step_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after the initial search
profile_params.step_growth_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after successfully extending the profile

% set bounds for optimizing when profiling the other parameters
profile_params.lb = lb;
profile_params.ub = ub;

% specify parameter ranges for bounds on profiling
profile_params.para_ranges = [profile_params.lb,profile_params.ub];

profile_params.opts = optim_opts;

%% objfn_constants
objfn_constants.fn = fn;
objfn_constants.fn_opts = fn_opts;
objfn_constants.weights = [1;1;1];

%% perform profile
profiles = performProfile(files,objfn_constants,profile_params,options);

%% save
save("data/" + file_name,"profiles")

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../ODEFittingFns/")
rmpath("../ODEFitting/")
