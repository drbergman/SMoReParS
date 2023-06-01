function P = optimizeSMParsFromABM(files,p0,fn,fn_opts,lb,ub,optim_opts,weights,input_opts)

opts = defaultOptimizeSMParsFromABMOptions;
if nargin >= 9 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

load(files.data_file,"t","D","C");
cohort_size = size(D,2:ndims(D)); % size of D is [n conditions, cohort_size]
D = reshape(D,size(D,1),[]); % string out all the cohorts along the 2nd dim
n_abm_vecs = size(D,2); % number of ABM parameter vectors used
m = size(D,1); % number of conditions used

assert(isempty(C) || (numel(C)==m && numel(weights)==m)) % make sure that these match if there are any conditions

npars = numel(p0);
P = zeros(npars,n_abm_vecs);
sz = [opts.n_starts,n_abm_vecs];
n_total = prod(sz);
if isfield(files,"previous_optim_file")
    load(files.previous_optim_file,"p_all","f_all");
    p_all = reshape(p_all,[npars,n_total]);
    f_all = reshape(f_all,sz);
    ind_to_run = find(f_all(:)==0);
else
    p_all = zeros([npars,n_total]);
    f_all = zeros(sz);
    ind_to_run = 1:n_total;
end
num_to_run = numel(ind_to_run);

last_save_time = tic;

if ~opts.force_serial
    if isempty(gcp('nocreate'))
        ppool = parpool("Processes");
    elseif isa(gcp,"parallel.ThreadPool")
        delete(gcp('nocreate'))
        ppool = parpool("Processes");
    else
        ppool = gcp;
    end

    %% This submits them all using parfeval
    FF(1:num_to_run) = parallel.FevalFuture;
    for i = 1:num_to_run
        [abm_par_ind,sample_ind] = ind2sub(sz,ind_to_run(i));
        if m==1
            d = D(abm_par_ind);
            F = @(p) rawError(p,t,d,fn,C{1},fn_opts,opts.raw_error_opts);
        else
            d = D(:,abm_par_ind);
            F = @(p) arrayfun(@(j) rawError(p,t,d(j),fn,C{j},fn_opts,opts.raw_error_opts),1:m)*weights;
        end
        if sample_ind>1
            p0_si = min(ub,max(lb,p0.*exp(.5*randn(npars,1))));
        else
            p0_si = p0;
        end
        FF(i) = parfeval(ppool,@() fmincon(F,p0_si,[],[],[],[],lb,ub,[],optim_opts),2);
    end

    for i = 1:num_to_run
        [temp_idx,p_temp,f_temp] = fetchNext(FF);
        idx = ind_to_run(temp_idx);
        p_all(:,idx) = p_temp;
        f_all(idx) = f_temp;

        if mod(i,round(.01*num_to_run))==0
            fprintf("Finished %d of %d.\n",i,num_to_run)
        end
        if mod(i,opts.save_every_iter)==0 && toc(last_save_time) > opts.save_every_sec
            save(opts.temp_profile_name,"p_all","f_all")
            fprintf("---------------Saved at iteration i = %d---------------\n",i)
            last_save_time = tic;
        end
    end

else
    p_all = reshape(p_all,[npars,sz]);
    for i = 1:num_to_run
        [abm_par_ind,sample_ind] = ind2sub(sz,ind_to_run(i));

        if sample_ind>1
            p0_si = min(ub,max(lb,p0.*exp(.5*randn(npars,1))));
        else
            p0_si = p0;
        end
        if m==1
            [p_all(:,sample_ind,abm_par_ind),f_all(sample_ind,abm_par_ind)] = fmincon(@(p) rawError(p,t,D(abm_par_ind),fn,C{1},fn_opts,opts.raw_error_opts),p0_si,[],[],[],[],lb,ub,[],optim_opts);
        else
            [p_all(:,sample_ind,abm_par_ind),f_all(sample_ind,abm_par_ind)] = fmincon(@(p) arrayfun(@(j) rawError(p,t,D(j,abm_par_ind),fn,C{j},fn_opts,opts.raw_error_opts),1:m)*weights,p0_si,[],[],[],[],lb,ub,[],optim_opts);
        end

        if mod(i,round(.01*num_to_run))==0
            fprintf("Finished %d of %d.\n",i,num_to_run)
        end
        if mod(i,opts.save_every_iter)==0 && toc(last_save_time) > opts.save_every_sec
            save(opts.temp_profile_name,"p_all","f_all")
            fprintf("---------------Saved at iteration i = %d---------------\n",i)
            last_save_time = tic;
        end
    end


end

p_all = reshape(p_all,[npars,sz]);
[~,best_case_inds] = min(f_all,[],1);
for i = 1:n_abm_vecs
    P(:,i) = p_all(:,best_case_inds(i),i);
end

P = reshape(P,[npars,cohort_size]);

end

function default_options = defaultOptimizeSMParsFromABMOptions

default_options.force_serial = false;
default_options.n_starts = 1; % number of times to start optimization, if >1 pick random perturbations of p0
default_options.save_every_iter = Inf; % how often (based on iterations) to save profile throughout (protects against errors and workers crashing)
default_options.save_every_sec = Inf; % how often (based on seconds passed) to save profile throughout (protects against errors and workers crashing)

temp_profile_name_format_spec = "data/temp_profile_%02d";
num = next_version_number(temp_profile_name_format_spec);
default_options.temp_profile_name = sprintf(temp_profile_name_format_spec,num); % how often to save profile throughout (protects against errors and workers crashing)

default_options.raw_error_opts = [];

end
