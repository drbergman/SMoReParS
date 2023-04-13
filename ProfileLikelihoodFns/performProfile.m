function out = performProfile(par_file,data_file,objfn_constants,profile_params,save_all_pars)


if ~isfield("objfn_constants","p_setup_fn")
    objfn_constants.p_setup_fn = @(p) p;
end

if nargin<5
    save_all_pars = false;
end

load(par_file,"P")
load(data_file,"t","D","C","cohort_size");
D = reshape(D,size(D,1),[]); % string out all the cohorts along the 2nd dim
n_abm_vecs = size(D,2); % number of ABM parameter vectors used
m = size(D,1); % number of conditions used

P = reshape(P,size(P,1),[]);
npars = size(P,1);

FF(1:n_abm_vecs) = parallel.FevalFuture;
t_start = tic;
out = cell(npars,n_abm_vecs);
for i = 1:n_abm_vecs
    d = D(:,i);
    p = P(:,i);
    % FF(i) = parfeval(@() profileLikelihood(p,t,d,C,objfn_constants,profile_params,save_all_pars),1);
    out(:,i) = profileLikelihood(p,t,d,C,objfn_constants,profile_params,save_all_pars);
end
fprintf("FevalQueue finished.\n")
%%
for i = 1:n_abm_vecs

    [idx,temp] = fetchNext(FF);
    out(:,idx) = temp;

    if mod(i,ceil(0.01*n_abm_vecs))==0
        temp = toc(t_start);
        fprintf("Finished %d after %s. ETR: %s\n",i,duration(0,0,temp),duration(0,0,temp/i * (n_abm_vecs-i)))
    end

end

out = reshape(out,[npars,cohort_size]);