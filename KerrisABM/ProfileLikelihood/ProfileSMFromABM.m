% This script sets up and calls the profile likelihood method for the ODE
% SM on a grid of ABM parameters. compare_every allows for comparing with a
% subset of ABM time points. This could be improved slightly by
% interpolating the ABM output at those times rather than the rounding I do
% to get tt. The function profileLikelihood computes the profile likelihood
% of all ODE model parameters at a given ABM parameter vector.

clearvars;
addpath("../../SurrogateModelFns/")
addpath("../../ProfileLikelihoodFns/")
addpath("~/Documents/MATLAB/myfunctions/")
addpath("../ODEFitting/")

opts.force_serial = true;
opts.profile_likelihood_options.save_all_pars = true;
opts.profile_likelihood_options.raw_error_opts.resample = true;
opts.profile_likelihood_options.raw_error_opts.t = 15:15:75;

%% objfn_constants
objfn_constants.fn = @computeTimeSeries;
objfn_constants.fn_opts.model_type = "logistic";
objfn_constants.weights = 1;

%% files
files.par_file = sprintf("../ODEFitting/data/OptimalParameters_%s.mat",objfn_constants.fn_opts.model_type);
files.data_file = "../PostAnalysis/data/summary.mat";
% files.previous_profile_file = "temp_profile.mat";

%% optimization opts
profile_params.opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%% setup SM parameters, bounds, etc
p = basePars(objfn_constants.fn_opts.model_type);
n_sm_pars = numel(p);

%% setup profile params
switch objfn_constants.fn_opts.model_type
    case "logistic"
        profile_params.initial_step_prop = .01*ones(n_sm_pars,1);
        profile_params.min_num_steps = 10*ones(n_sm_pars,1);
        profile_params.smallest_par_step = [1e-1;50]; % do not let the step size go below this as it steps towards the boundary/threshold
        profile_params.shrinking_factor = 0.9; % factor by which to shrink dx as it gets close to lower boundary
        profile_params.threshold = chi2inv(0.95,n_sm_pars); % compute threshold value for the parameter confidence intervals

        profile_params.secondary_step_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after the initial search
        profile_params.step_growth_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after successfully extending the profile

        % set bounds for optimizing when profiling the other parameters
        profile_params.lb = [0;0];
        profile_params.ub = [1;1e6];

        % specify parameter ranges for bounds on profiling
        profile_params.para_ranges = [0,1;     % r
            0,1e6];  % K

    case "von_bertalanffy"
        profile_params.initial_step_prop = .01*ones(n_sm_pars,1);
        profile_params.min_num_steps = 10*ones(n_sm_pars,1);
        profile_params.smallest_par_step = [1e-1;1e-1;1e-1]; % do not let the step size go below this as it steps towards the boundary/threshold
        profile_params.shrinking_factor = 0.9; % factor by which to shrink dx as it gets close to lower boundary
        profile_params.threshold = chi2inv(0.95,n_sm_pars); % compute threshold value for the parameter confidence intervals

        profile_params.secondary_step_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after the initial search
        profile_params.step_growth_factor = 2*ones(n_sm_pars,1); % factor by which to increase the step size after successfully extending the profile

        % set bounds for optimizing when profiling the other parameters
        profile_params.lb = [0;1;0];
        profile_params.ub = [100;1e3;100];
        profile_params.A = [-1,0,1]; % -alpha + beta <= 0 <==> alpha >= beta ==> population grows for small enough x
        profile_params.b = 0;

        % specify parameter ranges for bounds on profiling
        profile_params.para_ranges = [0,100;     % alpha
            1,1e3;  % nu
            0,100]; % beta for chemo activating apoptosis

    otherwise
        error("%s is an unspecified SM model.\n",objfn_constants.fn_opts.model_type);
end

%% perform profile
out = performProfile(files,objfn_constants,profile_params,opts);

save(sprintf("data/ProfileLikelihoods_%s.mat",objfn_constants.fn_opts.model_type),"out")

%% reset path
rmpath("../ODEFitting/")
rmpath("../../SurrogateModelFns/")
rmpath("../../ProfileLikelihoodFns/")


