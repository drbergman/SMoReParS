function P = optimizeSMParsFromABM(data_file,p0,fn,fn_opts,lb,ub,optim_opts,weights,input_opts)

opts = defaultOptimizeSMParsFromABMOptions;
if nargin >= 9 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

load(data_file,"t","D","C");
cohort_size = size(D,2:ndims(D)); % size of D is [n conditions, cohort_size]
D = reshape(D,size(D,1),[]); % string out all the cohorts along the 2nd dim
n_abm_vecs = size(D,2); % number of ABM parameter vectors used
m = size(D,1); % number of conditions used

assert(isempty(C) || (numel(C)==m && numel(weights)==m)) % make sure that these match if there are any conditions

npars = numel(p0);
P = zeros(npars,n_abm_vecs);

if ~opts.force_serial
    if isempty(gcp('nocreate'))
        ppool = parpool("Processes");
        % ppool = parpool("Threads");
    elseif isa(gcp,"parallel.ThreadPool")
    % elseif isa(gcp,"parallel.ProcessPool")
        delete(gcp('nocreate'))
        ppool = parpool("Processes");
        % ppool = parpool("Threads");
    else
        ppool = gcp;
    end

    %% This submits them all using parfeval
    FF(1:n_abm_vecs) = parallel.FevalFuture;

    ticBytes(ppool)
    for i = 1:n_abm_vecs
        if m==1
            d = D(i);
            F = @(p) rawError(p,t,d,fn,C{1},fn_opts,opts.raw_error_opts);
        else
            d = D(:,i);
            F = @(p) arrayfun(@(j) rawError(p,t,d(j),fn,C{j},fn_opts,opts.raw_error_opts),1:m)*weights;
        end
        FF(i) = parfeval(ppool,@() fmincon(F,p0,[],[],[],[],lb,ub,[],optim_opts),1);
    end

    for i = 1:n_abm_vecs
        [idx,temp] = fetchNext(FF);
        P(:,idx) = temp;

        if mod(i,round(.01*n_abm_vecs))==0
            fprintf("Finished %d of %d.\n",i,n_abm_vecs)
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
    for i = 1:n_abm_vecs
        if m==1
            P(:,i) = fmincon(@(p) rawError(p,t,D(i),fn,C{1},fn_opts,opts.raw_error_opts),p0,[],[],[],[],lb,ub,[],optim_opts);
        else
            % P(:,i) = arrayfun(@(j) fmincon(@(p) rawError(p,t,D(j,i),fn,C{j},fn_opts,opts.raw_error_opts),p0,[],[],[],[],lb,ub,[],optim_opts),1:m)*weights;
            P(:,i) = fmincon(@(p) arrayfun(@(j) rawError(p,t,D(j,i),fn,C{j},fn_opts,opts.raw_error_opts),1:m)*weights,p0,[],[],[],[],lb,ub,[],optim_opts);
        end

        if mod(i,round(.01*n_abm_vecs))==0
            fprintf("Finished %d of %d.\n",i,n_abm_vecs)
        end
        
    end
end

P = reshape(P,[npars,cohort_size]);

end

function default_options = defaultOptimizeSMParsFromABMOptions

default_options.force_serial = false;
default_options.raw_error_opts = [];

end
