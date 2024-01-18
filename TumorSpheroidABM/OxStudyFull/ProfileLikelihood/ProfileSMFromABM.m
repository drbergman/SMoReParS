
% This script sets up and calls the profile likelihood method for the ODE
% SM on a grid of ABM parameters. compare_every allows for comparing with a
% subset of ABM time points. This could be improved slightly by
% interpolating the ABM output at those times rather than the rounding I do
% to get tt. The function profileLikelihood computes the profile likelihood
% of all ODE model parameters at a given ABM parameter vector.

clearvars;
addpath("../../../SurrogateModelFns/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

addpath("~/Documents/MATLAB/myfunctions/")

file_name = "Profiles_SMFromABM_LMS_bounded";

files.optimal_parameters = "../ODEFitting/data/SMFitToABM_LMS_bounded.mat";

load(files.optimal_parameters,"cohort_name")

files.data = sprintf("../../data/%s/summary.mat",cohort_name);
% files.previous_profile_file = "data/temp_profile.mat";
% files.previous_profile_file = "ProfileLikelihoods.mat";

save_all_pars = true;
force_serial = false;
temp_profile_name = "data/temp_profile";
save_every_iter = 10; % wait at least this many iterations between saves
save_every_sec = 5*60; % wait at least this many seconds between saves

load("../ODEFitting/data/SMFitToData_LMS_bounded.mat","fixed_pars","lb","ub","model_type","optim_opts")
load("../ODEFitting/data/SMFitToData_LMS_bounded.mat","fn","fn_opts","sm")
if ~exist("sm","var")
    sm.fn = fn;
    sm.opts = fn_opts;
end
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
profile_params.weights = [1;1;1];

%% perform profile
profiles = performProfile(files,sm,profile_params,save_all_pars=save_all_pars,...
    force_serial=force_serial,temp_profile_name=temp_profile_name,...
    save_every_iter=save_every_iter,save_every_sec=save_every_sec);

%% save
save("data/" + file_name,"profiles")

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../SurrogateModelFns/")
rmpath("../ODEFitting/")
