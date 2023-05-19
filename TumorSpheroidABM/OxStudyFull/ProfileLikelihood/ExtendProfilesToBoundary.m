% This script is a one-off to extend profiles to the boundary of the ranges
% allowed for each parameter.

clearvars;

addpath("../../../ODEFittingFns/")
cohort_name = "cohort_2303301105";

% data_file = "../ODEFitting/data/ExperimentalData.mat";
data_file = sprintf("../../data/%s/summary_new.mat",cohort_name);
% profile_file = "data/Profiles_SMFromData";
profile_file = "data/Profiles_SMFromABM";

%% load profile
load(profile_file,"out")
load(data_file,"t","D","C","cohort_size");

%% set parameter ranges
n_sm_pars = size(out,1);
out = reshape(out,n_sm_pars,[]);
para_ranges = [0,20;     % lambda
               0,200;  % alpha
               0,1e4;      % K
               0,20;   % d in G1/S
               0,20;   % d in G2/M
               0,40]; % ec50 for chemo activating apoptosis
smallest_par_step = 1e-1*ones(n_sm_pars,1);
smallest_par_step(3) = 10; % use a larger min step for K
lb = [0;0;200;0;0;0];
ub = [20;200;1e4;20;20;40];
opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
D = reshape(D,size(D,1),[]); % string out all the cohorts along the 2nd dim

%% objfn_constants
fixed_hill_coefficient = 3;

objfn_constants.fn = @computeTimeSeries;
objfn_constants.fn_opts.phase_dependent_death = true; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?
objfn_constants.fn_opts.link_phase_death_rates = false; % whether to link the two phases death rates
objfn_constants.fn_opts.hill_activation = true; % if unlinked, use hill activation?

objfn_constants.weights = [1;1;1];
objfn_constants.p_setup_fn = @(p) [p(1:5);fixed_hill_coefficient;p(6)];

m = 3;

for i = 1:n_sm_pars
    for j = 1:size(out,2)
        F = @(p) arrayfun(@(k) rawError(objfn_constants.p_setup_fn(p),t,...
            D(k,j),objfn_constants.fn,C{k},objfn_constants.fn_opts),1:m)*objfn_constants.weights;
        if out{i,j}(i,1) < para_ranges(i,1) + smallest_par_step(i)
            x0 = out{i,j}(1:end-1,1);
            x0(i) = para_ranges(i,1);
            lb_temp = lb;
            ub_temp = ub;
            lb_temp(i) = x0(i);
            ub_temp(i) = x0(i);
            [x0_new,temp_val] = fmincon(F,x0,[],[],[],[],lb_temp,ub_temp,[],opts);
            out{i,j} = [[x0_new;temp_val],out{i,j}];
        end
        if out{i,j}(i,end) > para_ranges(i,2) - smallest_par_step(i)
            x0 = out{i,j}(1:end-1,end);
            x0(i) = para_ranges(i,2);
            lb_temp = lb;
            ub_temp = ub;
            lb_temp(i) = x0(i);
            ub_temp(i) = x0(i);
            [x0_new,temp_val] = fmincon(F,x0,[],[],[],[],lb_temp,ub_temp,[],opts);
            out{i,j} = [out{i,j},[x0_new;temp_val]];
        end
    end
end

%% save new profile
save(profile_file + "_extended","out")

%% reset path
rmpath("../../../ODEFittingFns/")
