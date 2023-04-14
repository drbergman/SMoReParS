function out = performProfile(files,objfn_constants,profile_params,save_all_pars,force_serial)


if ~isfield("objfn_constants","p_setup_fn")
    objfn_constants.p_setup_fn = @(p) p;
end

if nargin<4
    save_all_pars = false;
end

if nargin<5
    force_serial = false;
end

if ~isfield(profile_params,"A")
    profile_params.A = [];
    profile_params.b = [];
end    

load(files.par_file,"P")
load(files.data_file,"t","D","C","cohort_size");
D = reshape(D,size(D,1),[]); % string out all the cohorts along the 2nd dim
n_abm_vecs = size(D,2); % number of ABM parameter vectors used
m = size(D,1); % number of conditions used

P = reshape(P,size(P,1),[]);
npars = size(P,1);
t_start = tic;
if isfield(files,"previous_profile_file")
    load(files.previous_profile_file,"out")
    out = reshape(out,npars,n_abm_vecs);
else
    out = cell(npars,n_abm_vecs);
end

if ~force_serial
    FF(1:n_abm_vecs) = parallel.FevalFuture;

    if isempty(gcp('nocreate'))
        ppool = parpool("Processes");
    else
        ppool = gcp;
    end
    for i = 1:n_abm_vecs
        if any(cellfun(@isempty,out(:,i))) % then this one wasn't done (or somehow wasn't finished)
            d = D(:,i);
            p = P(:,i);
            FF(i) = parfeval(ppool,@() profileLikelihood(p,t,d,C,objfn_constants,profile_params,save_all_pars),1);
        else
            FF(i) = parfeval(@() "done",1); % if not empty, then just leave it be
        end
        % out(:,i) = profileLikelihood(p,t,d,C,objfn_constants,profile_params,save_all_pars);
    end
    fprintf("FevalQueue finished.\n")
    %%
    for i = 1:n_abm_vecs

        [idx,temp] = fetchNext(FF);
        if ~isequal(temp,"done")
            out(:,idx) = temp;
        end

        if mod(i,ceil(0.01*n_abm_vecs))==0
            temp = toc(t_start);
            fprintf("Finished %d after %s. ETR: %s\n",i,duration(0,0,temp),duration(0,0,temp/i * (n_abm_vecs-i)))
        end

    end
else
    for i = 1:n_abm_vecs
        if any(cellfun(@isempty,out(:,i))) % then this one wasn't done (or somehow wasn't finished)
            d = D(:,i);
            p = P(:,i);
            out(:,i) = profileLikelihood(p,t,d,C,objfn_constants,profile_params,save_all_pars);
        end
        if mod(i,ceil(0.01*n_abm_vecs))==0
            temp = toc(t_start);
            fprintf("Finished %d after %s. ETR: %s\n",i,duration(0,0,temp),duration(0,0,temp/i * (n_abm_vecs-i)))
        end
    end
end
out = reshape(out,[npars,cohort_size]);