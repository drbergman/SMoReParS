clearvars

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../ProfileLikelihoodFns/")
addpath("../../SurrogateModelFns/")
addpath("..")

files.optimal_parameters = "../SMFitting/data/OptimalParameters.mat";

n_sm_pars = 1;
force_serial = true;
save_all_pars = true;

%% setup profile params
profile_params.initial_step_prop = .01*ones(n_sm_pars,1);
profile_params.min_num_steps = 10*ones(n_sm_pars,1);
profile_params.shrinking_factor = 0.9; % factor by which to shrink dx as it gets close to lower boundary
profile_params.threshold = chi2inv(0.95,n_sm_pars); % compute threshold value for the parameter confidence intervals
profile_params.secondary_step_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after the initial search
profile_params.step_growth_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after successfully extending the profile
profile_params.smallest_par_step = 0.01; % do not let the step size go below this as it steps towards the boundary/threshold
profile_params.para_ranges = [-100,100];     % lambda

%% perform profile
profiles = performProfile(files,struct(),profile_params,force_serial=force_serial,...
    save_all_pars=save_all_pars);

mkdir data
save("data/Profiles.mat","profiles")



rmpath("../../ProfileLikelihoodFns/")
rmpath("../../SurrogateModelFns/")
rmpath("..")
