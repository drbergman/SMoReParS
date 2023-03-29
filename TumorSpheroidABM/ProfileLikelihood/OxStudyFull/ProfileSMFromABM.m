
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
cohort_name = "cohort_2303271138";

phase_dependent_death = false;

% load data and compute the standard error
if ~phase_dependent_death
    load("../../ODEFitting/OxStudyFull/data/OptimalParameters.mat","P")
else
    load("../../ODEFitting/OxStudyFull/data/OptimalParameters_phase_dependent_death.mat","P")
end
C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size","lattice_parameters");
chemo_dim = 8; % dimension along which chemo concentration varies; make sure this is still dim 8!!
if ~any(cohort_name==["cohort_2303231625","cohort_2303271138"])
    error("make sure the chemo concentration dim is still the 8th!")
end

Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"ode_state*");
load(sprintf("../../data/%s/output.mat",cohort_name),"ids");
load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");
t_abm = tracked.t;

%% sample from the ABM output
tt = t_abm;

%% setup Avg and Std so that they go [time,phase,death rate per uM,all other pars]
Avg = Sum.ode_state_average;
Avg = permute(Avg,[1,2,chemo_dim+2,setdiff(3:ndims(Avg),chemo_dim+2)]); % move chemo concentration dim to right after time and phase
Avg = reshape(Avg,numel(t_abm),2,size(Avg,3),[]); 

Std = Sum.ode_state_std;
Std = permute(Std,[1,2,chemo_dim+2,setdiff(3:ndims(Std),chemo_dim+2)]); % move chemo concentration dim to right after time and phase
Std = reshape(Std,numel(t_abm),2,size(Std,3),[]);

%%
threshold = chi2inv(0.95,3); % compute threshold value for the parameter confidence intervals

% specify parameter ranges for bounds on profiling
para_ranges = [0,10;     % lambda
               0,20;  % alpha
               0,5000;      % K
               0,2];   % d
min_par_step = [0.1;0.2;10;2e-2]; % d needs a minimum step because it can sometimes be very close to 0 at best fit
% set bounds for optimizing when profiling the other parameters
lb = [0;0;0;0];
ub = [Inf;Inf;1e4;2];
% opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
opts = optimset('Display','off');

% specify parameter names
para_names = {'\lambda', '\alpha', 'K', 'd'};
para_names_file_save = {'lambda', 'alpha', 'K', 'd'};

npars = size(P,1); % number of parameters

%%
n_abm_vecs = size(Avg,4);
FF(1:n_abm_vecs) = parallel.FevalFuture;
objfn_opts = struct("is_abm_data",true,"phase_dependent_death",phase_dependent_death);
chemo_vals = C.lattice_parameters(chemo_dim).values;
t_start = tic;
for i = 1:n_abm_vecs
% n_to_do = 8;
% time_taken = zeros(n_to_do,1);
% out = cell(4,n_to_do);
% for i = 1:n_to_do
    data = Avg(:,:,:,i);
    data_std = Std(:,:,:,i);
    data_std(data_std==0) = 1; % when the SD is zero, just compute the SM difference from the ABM rather than the z-score
    FF(i) = parfeval(@() profileLikelihood(P(:,i),t_abm,data,data_std,chemo_vals,phase_dependent_death,para_ranges,lb,ub,opts,threshold,min_par_step),1);
    % tic;
    % out(:,i) = profileLikelihood(P(:,i),t_abm,data,data_std,chemo_vals,phase_dependent_death,para_ranges,lb,ub,opts,threshold,min_par_step);
    % time_taken(i) = toc;

end

%%
out = cell(npars,n_abm_vecs);
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

