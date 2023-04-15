function P = optimizeSMParsFromABM(data_file,p,fn,fn_opts,lb,ub,opts,weights,force_serial)

if ~exist("force_serial","var")
    force_serial = false;
end

load(data_file,"t","D","C");
cohort_size = size(D,2:ndims(D)); % size of D is [n conditions, cohort_size]
D = reshape(D,size(D,1),[]); % string out all the cohorts along the 2nd dim
n_abm_vecs = size(D,2); % number of ABM parameter vectors used
m = size(D,1); % number of conditions used

assert(isempty(C) || (numel(C)==m && numel(weights)==m)) % make sure that these match if there are any conditions

npars = numel(p);
P = zeros(npars,n_abm_vecs);

if ~force_serial
    if isempty(gcp('nocreate'))
        ppool = parpool("Processes");
    else
        ppool = gcp;
    end
    FF(1:n_abm_vecs) = parallel.FevalFuture;

    for i = 1:n_abm_vecs
        if m==1
            F = @(p) rawError(p,t,D(i),fn,C{1},fn_opts);
        else
            d = D(:,i);
            F = @(p) arrayfun(@(j) rawError(p,t,d(j),fn,C{j},fn_opts),1:m)*weights;
        end
        FF(i) = parfeval(ppool,@() fmincon(F,p,[],[],[],[],lb,ub,[],opts),1);
    end

    for i = 1:n_abm_vecs
        [idx,temp] = fetchNext(FF);
        P(:,idx) = temp;
    end

else
    for i = 1:n_abm_vecs
        if m==1
            P(:,i) = fmincon(@(p) rawError(p,t,D(i),fn,C{1},fn_opts),p,[],[],[],[],lb,ub,[],opts);
        else
            P(:,i) = arrayfun(@(j) fmincon(@(p) rawError(p,t,D(j,i),fn,C{j},fn_opts),p,[],[],[],[],lb,ub,[],opts),1:m)*weights;
        end

        if mod(i,round(.01*n_abm_vecs))==0
            fprintf("Finished %d of %d.\n",i,n_abm_vecs)
        end
        
    end
end

P = reshape(P,[npars,cohort_size]);

