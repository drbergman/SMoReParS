
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
files.previous_profile_file = "ProfileLikelihoods.mat";

options.save_all_pars = true;
options.force_serial = false;
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

% %%
% C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size","lattice_parameters");
% chemo_dim = 8; % dimension along which chemo concentration varies; make sure this is still dim 8!!
% if ~any(cohort_name==["cohort_2303231625","cohort_2303271138"]) && ~isequal(C.lattice_parameters(chemo_dim).path,["chemo_pars","concentration"])
%     error("make sure the chemo concentration dim is still the 8th!")
% end
% 
% Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"count_*","state2_prop_*");
% load(sprintf("../../data/%s/output.mat",cohort_name),"ids","lattice_parameters");
% load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");
% tt = tracked.t';
% 
% %% setup profile params
% profile_params.initial_step_prop = .05*ones(6,1);
% profile_params.min_num_steps = 2*ones(6,1);
% profile_params.smallest_par_step = 1e-2*ones(6,1); % do not let the step size go below this as it steps towards the lower boundary
% profile_params.shrinking_factor = 0.9; % factor by which to shrink dx as it gets close to lower bound
% profile_params.threshold = chi2inv(0.95,6); % compute threshold value for the parameter confidence intervals
% 
% profile_params.min_par_step = [0.5;0.2;100;0.5;1;0.6]; % d needs a minimum step because it can sometimes be very close to 0 at best fit
% 
% % set bounds for optimizing when profiling the other parameters
% profile_params.lb = [0;0;0;0;0;0];
% profile_params.ub = [Inf;Inf;1e4;10;10;10];
% % opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
% profile_params.opts = optimset('Display','off');
% 
% % specify parameter ranges for bounds on profiling
% profile_params.para_ranges = [0,5;     % lambda
%                0,20;  % alpha
%                0,5000;      % K
%                0,5;   % d in G1/S
%                0,10;   % d in G2/M
%                0,6]; % ec50 for chemo activating apoptosis
% 
% %% setup Avg and Std so that they go [time,death rate per uM,all other pars]
% count = Sum.count_average;
% count = permute(count,[1,chemo_dim+1,setdiff(2:ndims(count),chemo_dim+1)]); % move chemo concentration dim to right after time
% count = reshape(count,numel(tt),size(count,2),[]); 
% 
% count_std = Sum.count_std;
% count_std = permute(count_std,[1,chemo_dim+1,setdiff(2:ndims(count_std),chemo_dim+1)]); % move chemo concentration dim to right after time
% count_std = reshape(count_std,numel(tt),size(count_std,2),[]);
% 
% state2_prop = Sum.state2_prop_mean;
% state2_prop = permute(state2_prop,[1,chemo_dim+1,setdiff(2:ndims(state2_prop),chemo_dim+1)]); % move chemo concentration dim to right after time
% state2_prop = reshape(state2_prop,numel(tt),size(state2_prop,2),[]); 
% 
% state2_prop_std = Sum.state2_prop_std;
% state2_prop_std = permute(state2_prop_std,[1,chemo_dim+1,setdiff(2:ndims(state2_prop_std),chemo_dim+1)]); % move chemo concentration dim to right after time
% state2_prop_std = reshape(state2_prop_std,numel(tt),size(state2_prop_std,2),[]);
% state2_prop_std(state2_prop_std==0) = 1; % if the STD is 0, just use the unnormalized difference
% 
% %%
% % specify parameter names
% para_names = {'\lambda', '\alpha', 'K', 'd'};
% para_names_file_save = {'lambda', 'alpha', 'K', 'd'};
% 
% 
% npars = size(P,1); % number of parameters
% 
% %%
% n_abm_vecs = size(count,3);
% FF(1:n_abm_vecs) = parallel.FevalFuture;
% objfn_constants.doses = C.lattice_parameters(chemo_dim).values;
% objfn_constants.tt = tt;
% objfn_constants.fn = @computeTimeSeries;
% t_start = tic;
% out = cell(npars,n_abm_vecs);
% save_all_pars = false;
% for i = 1:n_abm_vecs
%     vals = cat(3,count(:,:,i),state2_prop(:,:,i));
%     stds = cat(3,count_std(:,:,i),state2_prop_std(:,:,i));
%     p = P(:,i);
%     FF(i) = parfeval(@() profileLikelihood(p,vals,stds,objfn_constants,profile_params,save_all_pars),1);
%     % out(:,i) = profileLikelihood(P(:,i),vals,stds,objfn_constants,profile_params,save_all_pars);
% end
% fprintf("FevalQueue finished.\n")
% %%
% for i = 1:n_abm_vecs
% 
%     [idx,temp] = fetchNext(FF);
%     out(:,idx) = temp;
% 
%     if mod(i,100)==0
%         temp = toc(t_start);
%         fprintf("Finished %d after %s. ETR: %s\n",i,duration(0,0,temp),duration(0,0,temp/i * (n_abm_vecs-i)))
%     end
% 
% end

save("ProfileLikelihoods.mat","out")

rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../ODEFittingFns/")
rmpath("../ODEFitting/")

