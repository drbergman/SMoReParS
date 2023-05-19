
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

% folder to store plots and text files
cohort_name = "cohort_2303301105";

files.par_file = "../ODEFitting/data/OptimalParameters_UnLinkedHill.mat";
files.data_file = sprintf("../../data/%s/summary_new.mat",cohort_name);
% files.previous_profile_file = "data/temp_profile.mat";
% files.previous_profile_file = "ProfileLikelihoods.mat";

options.save_all_pars = true;
options.force_serial = true;
options.temp_profile_name = "data/temp_profile";
options.save_every_iter = 100; % wait at least this many iterations between saves
options.save_every_sec = 10*60; % wait at least this many seconds between saves

n_sm_pars = 6;

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
profile_params.lb = [0;0;200;0;0;0];
profile_params.ub = [20;200;1e4;20;20;40];

profile_params.opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

% specify parameter ranges for bounds on profiling
profile_params.para_ranges = [0,20;     % lambda
               0,200;  % alpha
               0,1e4;      % K
               0,20;   % d in G1/S
               0,20;   % d in G2/M
               0,40]; % ec50 for chemo activating apoptosis

%% objfn_constants
fixed_hill_coefficient = 3;

objfn_constants.fn = @computeTimeSeries;
objfn_constants.fn_opts.phase_dependent_death = true; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?
objfn_constants.fn_opts.link_phase_death_rates = false; % whether to link the two phases death rates
objfn_constants.fn_opts.hill_activation = true; % if unlinked, use hill activation?

objfn_constants.weights = [1;1;1];
objfn_constants.p_setup_fn = @(p) [p(1:5);fixed_hill_coefficient;p(6)];
%% perform profile
out = performProfile(files,objfn_constants,profile_params,options);

%% save
% save("data/Profiles_SMFromABM.mat","out")

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../ODEFittingFns/")
rmpath("../ODEFitting/")

