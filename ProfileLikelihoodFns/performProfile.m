function out = performProfile(files,objfn_constants,profile_params,input_opts)

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