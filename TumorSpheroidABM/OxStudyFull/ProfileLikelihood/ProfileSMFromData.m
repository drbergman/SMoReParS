
% This script sets up and calls the profile likelihood method for
% identifiability for the ODE pars from the data. it then plots the
% combinations

clearvars;
addpath("../../../ODEFittingFns/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

addpath("~/Documents/MATLAB/myfunctions/")

file_name = "Profiles_SMFromData_LMS_bounded";

files.optimal_parameters = "../ODEFitting/data/SMFitToData_LMS_bounded.mat";
files.data = "../ODEFitting/data/ExperimentalData.mat";
% files.previous_profile_file = "ProfileLikelihoods.mat";

load(files.optimal_parameters,"fixed_pars","fn","lb","ub","fn_opts","model_type","optim_opts")

options.profile_likelihood_options.save_all_pars = true;
options.force_serial = true; % no benefit to running this in parallel (only for doing this across multiple ABM parameter vectors)
% options.temp_profile_name = "data/temp_profile";
% options.save_every_iter = 100; % wait at least this many iterations between saves
% options.save_every_sec = 10*60; % wait at least this many seconds between saves

optim_opts.Display = "off";
[p,~,~,~] = fixParameters(model_type,fixed_pars);

%% process fixed pars
n_sm_pars = numel(p);

%% setup profile params
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
profile_params.para_ranges = [profile_params.lb,profile_params.ub];

profile_params.opts = optim_opts;

%% objfn_constants
objfn_constants.fn = fn;
objfn_constants.fn_opts = fn_opts;
objfn_constants.weights = [1;1;1];

%% perform profile
profiles = performProfile(files,objfn_constants,profile_params,options);

%% save the output
save("data/" + file_name,"profiles")

rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../ODEFittingFns/")
rmpath("../ODEFitting/")


