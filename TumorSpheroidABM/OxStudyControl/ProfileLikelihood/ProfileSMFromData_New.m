
% This script sets up and calls the profile likelihood method for
% identifiability for the ODE pars from the data. it then plots the
% combinations

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../../../ODEFittingFns/")
addpath("../ODEFitting/")

files.optimal_parameters = "../ODEFitting/data/SMFittoData_New.mat";
files.data = "../ODEFitting/data/ExperimentalData_New.mat";
% files.previous_profile_file = "ProfileLikelihoods.mat";

options.profile_likelihood_options.save_all_pars = true;
options.force_serial = true; % no benefit to running this in parallel (only for doing this across multiple ABM parameter vectors)

%% load data
load(files.optimal_parameters,"fn","lb","ub","fn_opts","optim_opts")
[p,~,~] = basePars();
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
objfn_constants.weights = 1;

%%
% specify parameter names
para_names = {'\lambda', '\alpha', 'K'};

%% compute profile likelihoods
profiles = performProfile(files,objfn_constants,profile_params,options);

%% save the output
save("data/Profiles_SMFromData_New.mat","profiles");
% ProfileLikelihoods_DataRestricted.mat used better bounds for the
% parameters in fitting. Use that.

%%  plot parameter combinations
figure;
ci=0;
for i = 1:n_sm_pars-1
    for j = i+1:n_sm_pars
        ci = ci+1;
        subplot(2,n_sm_pars,r2c(2,n_sm_pars,[1,ci]))
        plot(profiles{i}(i,:),profiles{i}(j,:))
        xlabel(para_names{i});ylabel(para_names{j})
        subplot(2,3,r2c(2,3,[2,ci]))
        plot(profiles{j}(j,:),profiles{j}(i,:))
        xlabel(para_names{j});ylabel(para_names{i})
    end
end

%%  plot parameter combinations
figure;
ci=0;
for i = 1:n_sm_pars-1
    for j = i+1:n_sm_pars
        ci = ci+1;
        subplot(1,n_sm_pars,r2c(1,n_sm_pars,[1,ci])); hold on
        plot(profiles{i}(i,:),profiles{i}(j,:)) % parameter j values as i was profiled
        xlabel(para_names{i});ylabel(para_names{j})
        plot(profiles{j}(i,:),profiles{j}(j,:)) % parameter i values as j was profiled (plotted as para i vs para j so it shares the same axes)
    end
end

%%
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../ODEFittingFns/")
rmpath("../ODEFitting/")
