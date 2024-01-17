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

cohort_name = "cohort_230124175743017";

files.par_file = "../ODEFitting/data/OptimalParameters_Using_optimizeSMParsFromABM.mat";
files.data_file = sprintf("../../data/%s/summary.mat",cohort_name);
% files.previous_profile_file = "temp_profile.mat";

save_all_pars = true;
force_serial = true;

n_sm_pars = 3;

%% setup profile params
profile_params.initial_step_prop = .01*ones(n_sm_pars,1);
profile_params.min_num_steps = 10*ones(n_sm_pars,1);
profile_params.smallest_par_step = [1e-1;1e-1;1e-1]; % do not let the step size go below this as it steps towards the boundary/threshold
profile_params.shrinking_factor = 0.9; % factor by which to shrink dx as it gets close to lower boundary
profile_params.threshold = chi2inv(0.95,n_sm_pars); % compute threshold value for the parameter confidence intervals

profile_params.secondary_step_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after the initial search
profile_params.step_growth_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after successfully extending the profile

% set bounds for optimizing when profiling the other parameters
profile_params.lb = [0;0;0];
profile_params.ub = [Inf;Inf;1e4];

profile_params.opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

% specify parameter ranges for bounds on profiling
profile_params.para_ranges = [0,100;     % lambda
               0,100;  % alpha
               0,1e4]; % K

%% objfn_constants
objfn_constants.fn = @computeTimeSeries;
objfn_constants.fn_opts = [];
objfn_constants.weights = 1;
out = performProfile(files,objfn_constants,profile_params,save_all_pars,force_serial);

save("data/ProfileLikelihoods.mat","out")

rmpath("../ODEFitting/")
rmpath("../../../SurrogateModelFns/")
rmpath("../../../ProfileLikelihoodFns/")


