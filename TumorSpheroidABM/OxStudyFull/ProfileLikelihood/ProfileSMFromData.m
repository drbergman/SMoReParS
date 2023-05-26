
% This script sets up and calls the profile likelihood method for
% identifiability for the ODE pars from the data. it then plots the
% combinations

clearvars;
addpath("../../../ODEFittingFns/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

addpath("~/Documents/MATLAB/myfunctions/")

files.par_file = "../ODEFitting/data/ODEFitToData.mat";
files.data_file = "../ODEFitting/data/ExperimentalData.mat";
% files.previous_profile_file = "ProfileLikelihoods.mat";

options.save_all_pars = true;
options.force_serial = true;
% options.temp_profile_name = "data/temp_profile";
% options.save_every_iter = 100; % wait at least this many iterations between saves
% options.save_every_sec = 10*60; % wait at least this many seconds between saves

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

%% save the output
save("data/Profiles_SMFromData.mat","out")

rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../ODEFittingFns/")
rmpath("../ODEFitting/")
% %%  plot parameter combinations
% figure;
% ci=0;
% for i = 1:N-1
%     for j = i+1:N
%         ci = ci+1;
%         subplot(2,N,r2c(2,N,[1,ci]))
%         plot(out{i}(i,:),out{i}(j,:))
%         xlabel(para_names{i});ylabel(para_names{j})
%         subplot(2,3,r2c(2,3,[2,ci]))
%         plot(out{j}(j,:),out{j}(i,:))
%         xlabel(para_names{j});ylabel(para_names{i})
%     end
% end
% 
% %%  plot parameter combinations
% figure;
% ci=0;
% for i = 1:N-1
%     for j = i+1:N
%         ci = ci+1;
%         subplot(1,N,r2c(1,N,[1,ci])); hold on
%         plot(out{i}(i,:),out{i}(j,:)) % parameter j values as i was profiled
%         xlabel(para_names{i});ylabel(para_names{j})
%         plot(out{j}(i,:),out{j}(j,:)) % parameter i values as j was profiled (plotted as para i vs para j so it shares the same axes)
%         xlabel(para_names{j});ylabel(para_names{i})
%     end
% end

%%

rmpath("../ODEFitting/")
rmpath("../../../ODEFittingFns/")
rmpath("../../../ProfileLikelihoodFns/")

