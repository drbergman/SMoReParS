% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;

make_save = false;
addpath("../../../ODEFittingFns/")
addpath("~/Documents/MATLAB/myfunctions/")

opts.force_serial = false;
opts.raw_error_opts.assume_independent_time_series = true; % assume that the two time series have diagonal covariance matrices at each time point

cohort_name = "cohort_2303301105";
data_file = sprintf("../../data/%s/summary_new.mat",cohort_name);

p = zeros(4,1);

p(1) = 24/19; % lambda
p(2) = 24/5; % alpha
p(3) = 1e3; % K
p(4) = 0.1; % chemo-induced death rate per uM of drug

p_unlinked = [1.9;1.9]; % chemo-induced death rates if they are unlinked
p_hill = [3;.1]; % [hill coefficient ; EC50] if unlinked

fn = @computeTimeSeries;
fn_opts.phase_dependent_death = true; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?
fn_opts.link_phase_death_rates = false; % whether to link the two phases death rates
fn_opts.hill_activation = true; % if unlinked, use hill activation?

ub = [Inf;Inf;1e4;10]; % running this once showed that beyond 0.5, delta just decreases populations too fast

optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

weights = [1;1;1];

%% finish setting up bounds
if ~fn_opts.link_phase_death_rates
    p = [p(1:3);p_unlinked];
    ub(5) = 10;
    if fn_opts.hill_activation
        p = [p;p_hill];
        ub(6:7) = [3;Inf];
    end
end

npars = length(p);
lb = zeros(npars,1);

if ~fn_opts.link_phase_death_rates
    if fn_opts.hill_activation
        lb(6) = 3;
    end
end

%% optimize sm pars
P = optimizeSMParsFromABM(data_file,p,fn,fn_opts,lb,ub,optim_opts,weights,opts);

if make_save
    save("data/OptimalParameters.mat","P")
end

rmpath("../../../ODEFittingFns/")

% %% load ABM data
% C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size","lattice_parameters");
% chemo_dim = 8; % dimension along which chemo concentration varies; make sure this is still dim 8!!
% if ~any(cohort_name==["cohort_2303231625","cohort_2303271138","cohort_2303301105"])
%     error("make sure the chemo concentration dim is still the 8th!")
% end
% 
% Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"count_*","state2_prop_*");
% load(sprintf("../../data/%s/output.mat",cohort_name),"ids");
% load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");
% t_abm = tracked.t;
% 
% %% set chemo val
% % chemo_val = zeros(size(ids));
% 
% 
% %% setup Avg and Std so that they go [time,death rate per uM,all other pars]
% count = Sum.count_average;
% count = permute(count,[1,chemo_dim+1,setdiff(2:ndims(count),chemo_dim+1)]); % move chemo concentration dim to right after time
% count = reshape(count,numel(t_abm),size(count,2),[]); 
% 
% count_std = Sum.count_std;
% count_std = permute(count_std,[1,chemo_dim+1,setdiff(2:ndims(count_std),chemo_dim+1)]); % move chemo concentration dim to right after time
% count_std = reshape(count_std,numel(t_abm),size(count_std,2),[]);
% 
% state2_prop = Sum.state2_prop_mean;
% state2_prop = permute(state2_prop,[1,chemo_dim+1,setdiff(2:ndims(state2_prop),chemo_dim+1)]); % move chemo concentration dim to right after time
% state2_prop = reshape(state2_prop,numel(t_abm),size(state2_prop,2),[]); 
% 
% state2_prop_std = Sum.state2_prop_std;
% state2_prop_std = permute(state2_prop_std,[1,chemo_dim+1,setdiff(2:ndims(state2_prop_std),chemo_dim+1)]); % move chemo concentration dim to right after time
% state2_prop_std = reshape(state2_prop_std,numel(t_abm),size(state2_prop_std,2),[]);
% 
% %%
% P = zeros([npars,size(count,3)]); % one parameter vector that works for each concentration
% 
% %% find the best fit pars; since the purpose of this is to give ranges for the ODE parameters for profiling later, the obj fun is not so critical
% nconc = C.cohort_size(chemo_dim);
% doses = C.lattice_parameters(chemo_dim).values;
% 
% if isempty(gcp('nocreate'))
%     ppool = parpool("Processes");
% else
%     ppool = gcp;
% end
% FF(1:size(P,2)) = parallel.FevalFuture;
% 
% state2_prop_std(state2_prop_std==0) = 1; % if the STD is 0, just use the unnormalized difference
% for i = 1:size(P,2)
%     temp_count = count(:,:,i);
%     temp_count_std = count_std(:,:,i);
%     temp_state2_prop = state2_prop(:,:,i);
%     temp_state2_prop_std = state2_prop_std(:,:,i);
%     F = @(p) arrayfun(@(j) rawError(p,t_abm,[temp_count(:,j),temp_state2_prop(:,j)],[temp_count_std(:,j),temp_state2_prop_std(:,j)],fn,doses(j),fn_opts),1:3)*weights;
%     FF(i) = parfeval(@() fmincon(F,p,[],[],[],[],lb,ub,[],optim_opts),1);
% end
% 
% for i = 1:size(P,2)
% 
%     [idx,temp] = fetchNext(FF);
%     P(:,idx) = temp;
% 
%     if mod(i,100)==0
%         fprintf("Finished %d.\n",i)
%     end
% 
% end