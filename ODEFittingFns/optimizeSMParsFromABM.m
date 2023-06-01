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
n_to_run = numel(ind_to_run);

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
    FF(1:n_to_run) = parallel.FevalFuture;
    ticBytes(ppool)
    for ii = 1:n_to_run
        [i,si] = ind2sub(sz,ind_to_run(ii));
        if m==1
            d = D(i);
            F = @(p) rawError(p,t,d,fn,C{1},fn_opts,opts.raw_error_opts);
        else
            d = D(:,i);
            F = @(p) arrayfun(@(j) rawError(p,t,d(j),fn,C{j},fn_opts,opts.raw_error_opts),1:m)*weights;
        end
        if si>1
            p0_si = min(ub,max(lb,p0.*exp(.5*randn(npars,1))));
        else
            p0_si = p0;
        end
        FF(ii) = parfeval(ppool,@() fmincon(F,p0_si,[],[],[],[],lb,ub,[],optim_opts),2);
    end

    for i = 1:n_to_run
        [temp_idx,p_temp,f_temp] = fetchNext(FF);
        idx = ind_to_run(temp_idx);
        p_all(:,idx) = p_temp;
        f_all(idx) = f_temp;

        if mod(i,round(.01*n_to_run))==0
            fprintf("Finished %d of %d.\n",i,n_to_run)
        end
        if mod(i,opts.save_every_iter)==0 && toc(last_save_time) > opts.save_every_sec
            save(opts.temp_profile_name,"p_all","f_all")
            fprintf("---------------Saved at iteration i = %d---------------\n",i)
            last_save_time = tic;
        end
    end
    tocBytes(ppool)

    %% this splits it all up into batches to help with sending too much data across the workerss at once and crashing matlab? I think that's why it crashes in the OxStudyFull case
    % max_in_batch = 400;
    % n_batches = ceil(n_abm_vecs/max_in_batch);
    % ticBytes(ppool)
    % 
    % for bi = 1:n_batches
    %     clear FF
    %     i_in_batch_start = 1+max_in_batch*(bi-1);
    %     i_in_batch_end = min(n_abm_vecs,max_in_batch*bi);
    %     i_in_batch = i_in_batch_start:i_in_batch_end;
    %     n_in_batch = length(i_in_batch);
    %     FF(1:n_in_batch) = parallel.FevalFuture;
    % 
    %     for i = 1:n_in_batch
    %         if m==1
    %             d = D(i_in_batch(i));
    %             F = @(p) rawError(p,t,d,fn,C{1},fn_opts,opts.raw_error_opts);
    %         else
    %             d = D(:,i_in_batch(i));
    %             F = @(p) arrayfun(@(j) rawError(p,t,d(j),fn,C{j},fn_opts,opts.raw_error_opts),1:m)*weights;
    %         end
    %         FF(i) = parfeval(ppool,@() fmincon(F,p0,[],[],[],[],lb,ub,[],optim_opts),1);
    %     end
    % 
    %     for i = 1:n_in_batch
    %         [idx,temp] = fetchNext(FF);
    %         P(:,i_in_batch(idx)) = temp;
    % 
    %         if mod(i_in_batch(i),round(.01*n_abm_vecs))==0
    %             fprintf("Finished %d of %d.\n",i_in_batch(i),n_abm_vecs)
    %         end
    %     end
    % end
    % tocBytes(ppool)

else
    p_all = reshape(p_all,[npars,sz]);
    for ii = 1:n_to_run
        [i,si] = ind2sub(sz,ind_to_run(ii));

        if si>1
            p0_si = min(ub,max(lb,p0.*exp(.5*randn(npars,1))));
        else
            p0_si = p0;
        end
        if m==1
            [p_all(:,si,i),f_all(si,i)] = fmincon(@(p) rawError(p,t,D(i),fn,C{1},fn_opts,opts.raw_error_opts),p0_si,[],[],[],[],lb,ub,[],optim_opts);
        else
            [p_all(:,si,i),f_all(si,i)] = fmincon(@(p) arrayfun(@(j) rawError(p,t,D(j,i),fn,C{j},fn_opts,opts.raw_error_opts),1:m)*weights,p0_si,[],[],[],[],lb,ub,[],optim_opts);
        end

        if mod(ii,round(.01*n_to_run))==0
            fprintf("Finished %d of %d.\n",ii,n_to_run)
        end
        if mod(ii,opts.save_every_iter)==0 && toc(last_save_time) > opts.save_every_sec
            save(opts.temp_profile_name,"p_all","f_all")
            fprintf("---------------Saved at iteration i = %d---------------\n",ii)
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
