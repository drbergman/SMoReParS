% This script sets up and calls the profile likelihood method for the ODE
% SM on a grid of ABM parameters. compare_every allows for comparing with a
% subset of ABM time points. This could be improved slightly by
% interpolating the ABM output at those times rather than the rounding I do
% to get tt. The function profileLikelihood computes the profile likelihood
% of all ODE model parameters at a given ABM parameter vector.

clearvars;
addpath("../../ODEFittingFns/")
addpath("../../ProfileLikelihoodFns/")
addpath("~/Documents/MATLAB/myfunctions/")
addpath("../ODEFitting/")

par_file = "../ODEFitting/data/OptimalParameters.mat";
data_file = "../PostAnalysis/summary.mat";

n_sm_pars = 3;

%% setup profile params
profile_params.initial_step_prop = .05*ones(n_sm_pars,1);
profile_params.min_num_steps = 2*ones(n_sm_pars,1);
profile_params.smallest_par_step = [1e-3;1e-3;1e-5]; % do not let the step size go below this as it steps towards the lower boundary
profile_params.shrinking_factor = 0.9; % factor by which to shrink dx as it gets close to lower bound
profile_params.threshold = chi2inv(0.95,n_sm_pars); % compute threshold value for the parameter confidence intervals

profile_params.min_par_step = [0.1;0.01;0.001]; % d needs a minimum step because it can sometimes be very close to 0 at best fit

% set bounds for optimizing when profiling the other parameters
profile_params.lb = [0;0;0];
profile_params.ub = [Inf;1;Inf];
profile_params.opts = optimset('Display','off');

% specify parameter ranges for bounds on profiling
profile_params.para_ranges = [0,100;     % alpha
               0,1;  % theta
               0,30]; % beta for chemo activating apoptosis

%% objfn_constants
objfn_constants.fn = @computeTimeSeries;
objfn_constants.fn_opts = [];
objfn_constants.weights = 1;
out = performProfile(par_file,data_file,objfn_constants,profile_params);

save("data/ProfileLikelihoods.mat","out")

rmpath("../ODEFitting/")
rmpath("../../ODEFittingFns/")
rmpath("../../ProfileLikelihoodFns/")


