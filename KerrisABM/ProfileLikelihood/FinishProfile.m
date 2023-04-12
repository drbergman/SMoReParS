
% This script continues the profiling started with ProfileSMFromABM.m in
% the event that that script was cut short (and half-finished output was
% saved).

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../ODEFitting/OxStudyFull/")

% folder to store plots and text files
cohort_name = "cohort_2303301105";
load("../../ODEFitting/OxStudyFull/data/OptimalParameters_UnLinkedHill.mat");

objfn_constants.hill_coefficient = 3;

fn = @computeTimeSeries;
objfn_constants.fn_opts.phase_dependent_death = true; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?
objfn_constants.fn_opts.link_phase_death_rates = false; % whether to link the two phases death rates
objfn_constants.fn_opts.hill_activation = true; % if unlinked, use hill activation?

objfn_constants.weights = [1;1;1];

C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size","lattice_parameters");
chemo_dim = 8; % dimension along which chemo concentration varies; make sure this is still dim 8!!
if ~any(cohort_name==["cohort_2303231625","cohort_2303271138"]) && ~isequal(C.lattice_parameters(chemo_dim).path,["chemo_pars","concentration"])
    error("make sure the chemo concentration dim is still the 8th!")
end

Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"count_*","state2_prop_*");
load(sprintf("../../data/%s/output.mat",cohort_name),"ids","lattice_parameters");
load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");
tt = tracked.t';

%% setup profile params
profile_params.initial_step_prop = .05*ones(6,1);
profile_params.min_num_steps = 2*ones(6,1);
profile_params.smallest_par_step = 1e-2*ones(6,1); % do not let the step size go below this as it steps towards the lower boundary
profile_params.shrinking_factor = 0.9; % factor by which to shrink dx as it gets close to lower bound
profile_params.threshold = chi2inv(0.95,6); % compute threshold value for the parameter confidence intervals

profile_params.min_par_step = [0.5;0.2;100;0.5;1;0.6]; % d needs a minimum step because it can sometimes be very close to 0 at best fit

% set bounds for optimizing when profiling the other parameters
profile_params.lb = [0;0;0;0;0;0];
profile_params.ub = [Inf;Inf;1e4;10;10;10];
% opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
profile_params.opts = optimset('Display','off');

% specify parameter ranges for bounds on profiling
profile_params.para_ranges = [0,5;     % lambda
               0,20;  % alpha
               0,5000;      % K
               0,5;   % d in G1/S
               0,10;   % d in G2/M
               0,6]; % ec50 for chemo activating apoptosis

%% setup Avg and Std so that they go [time,death rate per uM,all other pars]
count = Sum.count_average;
count = permute(count,[1,chemo_dim+1,setdiff(2:ndims(count),chemo_dim+1)]); % move chemo concentration dim to right after time
count = reshape(count,numel(tt),size(count,2),[]); 

count_std = Sum.count_std;
count_std = permute(count_std,[1,chemo_dim+1,setdiff(2:ndims(count_std),chemo_dim+1)]); % move chemo concentration dim to right after time
count_std = reshape(count_std,numel(tt),size(count_std,2),[]);

state2_prop = Sum.state2_prop_mean;
state2_prop = permute(state2_prop,[1,chemo_dim+1,setdiff(2:ndims(state2_prop),chemo_dim+1)]); % move chemo concentration dim to right after time
state2_prop = reshape(state2_prop,numel(tt),size(state2_prop,2),[]); 

state2_prop_std = Sum.state2_prop_std;
state2_prop_std = permute(state2_prop_std,[1,chemo_dim+1,setdiff(2:ndims(state2_prop_std),chemo_dim+1)]); % move chemo concentration dim to right after time
state2_prop_std = reshape(state2_prop_std,numel(tt),size(state2_prop_std,2),[]);
state2_prop_std(state2_prop_std==0) = 1; % if the STD is 0, just use the unnormalized difference

%%
% specify parameter names
para_names = {'\lambda', '\alpha', 'K', 'd'};
para_names_file_save = {'lambda', 'alpha', 'K', 'd'};


npars = size(P,1); % number of parameters

%%
n_abm_vecs = size(count,3);
FF(1:n_abm_vecs) = parallel.FevalFuture;
objfn_constants.doses = C.lattice_parameters(chemo_dim).values;
objfn_constants.tt = tt;
objfn_constants.fn = @computeTimeSeries;
t_start = tic;

load("temp_profile.mat","out")

save_all_pars = false;
for i = 1:n_abm_vecs
    if any(cellfun(@isempty,out(:,i))) % then this one wasn't done (or somehow wasn't finished)
        vals = cat(3,count(:,:,i),state2_prop(:,:,i));
        stds = cat(3,count_std(:,:,i),state2_prop_std(:,:,i));
        p = P(:,i);
        FF(i) = parfeval(@() profileLikelihood(p,vals,stds,objfn_constants,profile_params,save_all_pars),1);
    else
        FF(i) = parfeval(@() "done",1); % if not empty, then just leave it be
    end
    % out(:,i) = profileLikelihood(P(:,i),vals,stds,objfn_constants,profile_params,save_all_pars);
end
fprintf("FevalQueue finished.\n")
%%
for i = 1:n_abm_vecs

    [idx,temp] = fetchNext(FF);
    if ~isequal(temp,"done")
        out(:,idx) = temp;
    end

    if mod(i,100)==0
        temp = toc(t_start);
        fprintf("Finished %d after %s. ETR: %s\n",i,duration(0,0,temp),duration(0,0,temp/i * (n_abm_vecs-i)))
    end

end

% save("ProfileLikelihoods.mat","out")

rmpath("../../ODEFitting/OxStudyFull/")

