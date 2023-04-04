
% This script sets up and calls the profile likelihood method for the ODE
% SM on a grid of ABM parameters. compare_every allows for comparing with a
% subset of ABM time points. This could be improved slightly by
% interpolating the ABM output at those times rather than the rounding I do
% to get tt. The function profileLikelihood computes the profile likelihood
% of all ODE model parameters at a given ABM parameter vector.

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../ODEFitting/OxStudyFull/")

% folder to store plots and text files
cohort_name = "cohort_2303301105";
load("data/OptimalParameters_UnLinkedHill.mat");


fn = @computeTimeSeries;
objfn_data.fn_opts.phase_dependent_death = true; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?
objfn_data.fn_opts.link_phase_death_rates = false; % whether to link the two phases death rates
objfn_data.fn_opts.hill_activation = true; % if unlinked, use hill activation?

objfn_data.weights = [1;1;1];

C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size","lattice_parameters");
chemo_dim = 8; % dimension along which chemo concentration varies; make sure this is still dim 8!!
if ~any(cohort_name==["cohort_2303231625","cohort_2303271138"]) && ~isequal(lattice_parameters(chemo_dim).path,["chemo_pars","concentration"])
    error("make sure the chemo concentration dim is still the 8th!")
end

Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"count_*","state2_prop_*");
load(sprintf("../../data/%s/output.mat",cohort_name),"ids","lattice_parameters");
load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");
tt = tracked.t';

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
threshold = chi2inv(0.95,3); % compute threshold value for the parameter confidence intervals

% specify parameter ranges for bounds on profiling
profile_params.para_ranges = [0,5;     % lambda
               0,20;  % alpha
               0,5000;      % K
               0,5;   % d in G1/S
               0,10;   % d in G2/M
               0,6]; % ec50 for chemo activating apoptosis

profile_params.min_par_step = [0.05;0.02;10;0.05;0.1;0.06]; % d needs a minimum step because it can sometimes be very close to 0 at best fit
% set bounds for optimizing when profiling the other parameters
profile_params.lb = [0;0;0;0;0;0];
profile_params.ub = [Inf;Inf;1e4;10;10;10];
% opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
profile_params.opts = optimset('Display','off');

% specify parameter names
para_names = {'\lambda', '\alpha', 'K', 'd'};
para_names_file_save = {'lambda', 'alpha', 'K', 'd'};


npars = size(P,1); % number of parameters

%%
n_abm_vecs = size(Avg,4);
FF(1:n_abm_vecs) = parallel.FevalFuture;
objfn_data.doses = C.lattice_parameters(chemo_dim).values;
objfn_data.tt = tt;
objfn_data.fn = @computeTimeSeries;
t_start = tic;
out = cell(npars,n_abm_vecs);
save_all_pars = false;
for i = 1:n_abm_vecs
    objfn_data.count = count(:,:,i);
    objfn_data.count_std = count_std(:,:,i);
    objfn_data.state2_prop = state2_prop(:,:,i);
    objfn_data.state2_prop_std = state2_prop_std(:,:,i);
    % FF(i) = parfeval(@() profileLikelihood(P(:,i),t_abm,data,data_std,chemo_vals,phase_dependent_death,para_ranges,lb,ub,opts,threshold,min_par_step),1);
    out(:,i) = profileLikelihood(P(:,i),objfn_data,profile_params,save_all_pars);
end

%%
for i = 1:n_abm_vecs

    [idx,temp] = fetchNext(FF);
    out(:,idx) = temp;

    if mod(i,100)==0
        temp = toc(t_start);
        fprintf("Finished %d after %s. ETR: %s\n",i,duration(0,0,temp),duration(0,0,temp/i * (n_abm_vecs-i)))
    end

end

% save("ProfileLikelihoods.mat","out")

rmpath("~/Documents/MATLAB/myfunctions/")
addpath("../../ODEFitting/OxStudyFull/")

