function out = performProfile(files,objfn_constants,profile_params,input_opts)

% profiles each SM parameter for each "column" of the data array, Data 
% array being reshaped to a 2D array. 
% 
% files: struct with fields 
%   par_file: best fit SM parameters for each column of Data
%   data_file: the data in size [nconditions,cohort_size]
%   previous_profile_file (optional): partially completed profile
% objfn_constants: struct with fields that do vary by condition
%   p_setup_fn (optional): function that transforms SM parameter vector (can be used to fix a given SM parameter, e.g.)
%   fn: function used by rawError to generate SM output to compare against data
%   fn_opts: options to be used by fn for all conditions
%   weights: weights to use across conditions for objective function
% profile_params: struct with fields
%   lb: lower bound of SM parameters in calls to fmincon
%   ub: upper bound of SM parameters in calls to fmincon
%   A (optional): array to give A*x<=b constraint in fmincon
%   b (optional): vector to give A*x<=b constraint in fmincon
%   opts: options to pass in to fmincon
%   step_growth_factor (optional): factor to increase step size in Stage 2 of profiling
%   para_ranges: n_sm_pars x 2 array of extrema to be used in profiling for each parameter, i.e. profiles never force this parameter beyond this value
%   shrinking_factor: factor by which to shrink step sizes when approaching bounds in para_ranges
%   initial_step_prop: proportion of best parameter to move in Stage 1
%   smallest_par_step: minimum allowable step size for each parameter; Stage 2 halts when the step size is below this and the threshold is exceeded
%   min_num_steps: number of steps for each SM to take in Stage 1
%   threshold: chi2inv value for profile
%   secondary_step_factor: factor to increase an SM parameter step size when beginning Stage 2
%   step_growth_factor: factor to increase step size in Stage 2 after successfully extending profile
% input_opts (optional): struct with fields
%   force_serial=true: logical to force serial computation of profiles
%   save_all_pars=true: logical to control whether to save all parameter values (true) or only the current profile parameter and goodness-of-fit value (false)
%   save_every_iter=Inf: how often (based on iterations) to save profile throughout (protects against errors and workers crashing)
%   save_every_sec=Inf: how often (based on seconds passed) to save profile throughout (protects against errors and workers crashing)
%   temp_profile_name: how often to save profile throughout (protects against errors and workers crashing)

opts = defaultPerformProfileOptions;
if nargin==4 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

if ~isfield(objfn_constants,"p_setup_fn")
    objfn_constants.p_setup_fn = @(p) p;
end

if ~isfield(profile_params,"A")
    profile_params.A = [];
    profile_params.b = [];
end    

load(files.par_file,"P")
load(files.data_file,"t","D","C","cohort_size");
D = reshape(D,size(D,1),[]); % string out all the cohorts along the 2nd dim
n_abm_vecs = size(D,2); % number of ABM parameter vectors used

P = reshape(P,size(P,1),[]);
npars = size(P,1);
t_start = tic;
if isfield(files,"previous_profile_file")
    load(files.previous_profile_file,"out")
    out = reshape(out,npars,n_abm_vecs);
    ind_to_run = find(any(cellfun(@isempty,out),1));
else
    out = cell(npars,n_abm_vecs);
    ind_to_run = 1:n_abm_vecs;
end

num_to_run = length(ind_to_run);
fprintf("Going to run (# SM Parameters) x (# ABM Parameters) = %d x %d = %d profiles.\n",npars,num_to_run,npars*num_to_run);

if ~isfield(profile_params,"step_growth_factor")
    profile_params.step_growth_factor = ones(npars,1);
end

last_save_time = tic;

if ~opts.force_serial 
    %% run in parallel
    FF(1:num_to_run) = parallel.FevalFuture;

    if isempty(gcp('nocreate'))
        ppool = parpool("Processes");
    else
        ppool = gcp;
    end
    for i = 1:num_to_run
        current_ind = ind_to_run(i);
        d = D(:,current_ind);
        p = P(:,current_ind);
        FF(i) = parfeval(ppool,@() profileLikelihood(p,t,d,C,objfn_constants,profile_params,opts.save_all_pars),1);
    end
    fprintf("FevalQueue finished.\n")
    %% fetch profiles performed in parallel
    for i = 1:num_to_run

        [idx,temp] = fetchNext(FF);
        out(:,ind_to_run(idx)) = temp;

        if mod(i,ceil(0.001*num_to_run))==0
            temp = toc(t_start);
            fprintf("Finished %d of %d after %s. ETR: %s\n",i,num_to_run,duration(0,0,temp),duration(0,0,temp/i * (num_to_run-i)))
        end
        if mod(i,opts.save_every_iter)==0 && toc(last_save_time) > opts.save_every_sec
            save(opts.temp_profile_name,"out")
            fprintf("---------------Saved at iteration i = %d---------------\n",i)
            last_save_time = tic;
        end
    end
else
    %% run in serial
    for i = 1:num_to_run
        current_ind = ind_to_run(i);
        d = D(:,current_ind);
        p = P(:,current_ind);
        out(:,current_ind) = profileLikelihood(p,t,d,C,objfn_constants,profile_params,opts.save_all_pars);
        if mod(i,ceil(0.001*num_to_run))==0
            temp = toc(t_start);    
            fprintf("Finished %d of %d after %s. ETR: %s\n",i,num_to_run,duration(0,0,temp),duration(0,0,temp/i * (num_to_run-i)))
        end
        if mod(i,opts.save_every_iter)==0 && toc(last_save_time) > opts.save_every_sec
            save(opts.temp_profile_name,"out")
            fprintf("---------------Saved at iteration i = %d---------------\n",i)
            last_save_time = tic;
        end
    end
end
out = reshape(out,[npars,cohort_size]);

end

function default_options = defaultPerformProfileOptions

default_options.force_serial = true;
default_options.save_all_pars = true;
default_options.save_every_iter = Inf; % how often (based on iterations) to save profile throughout (protects against errors and workers crashing)
default_options.save_every_sec = Inf; % how often (based on seconds passed) to save profile throughout (protects against errors and workers crashing)

temp_profile_name_format_spec = "data/temp_profile_%02d";
num = next_version_number(temp_profile_name_format_spec);
default_options.temp_profile_name = sprintf(temp_profile_name_format_spec,num); % how often to save profile throughout (protects against errors and workers crashing)


end