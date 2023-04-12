function out = performProfile(par_file,data_file,objfn_constants)


if ~isfield("objfn_constants","p_setup_fn")
    objfn_constants.p_setup_fn = @(p) p;
end

load(par_file,"P")
load(data_file,"t","D","S","C","cohort_size","nsamps_per_parameter_vector","n_time_series","n_conditions");

if ~iscell(D)
    error("Expect the data to be a cell array where rows correspond to conditions and higher dimensions correspond to cohorts (e.g. if profiling against an ABM these would be all the ABM parameter vectors)")
end

P = reshape(P,size(P,1),[]);
npars = size(P,1);
n_abm_vecs = size(P,2);

FF(1:n_abm_vecs) = parallel.FevalFuture;
t_start = tic;
out = cell(npars,n_abm_vecs);
save_all_pars = false;
for i = 1:n_abm_vecs
    d = D(:,i);
    s = S(:,i);
    p = P(:,i);
    FF(i) = parfeval(@() profileLikelihood(p,t,d,s,objfn_constants,profile_params,save_all_pars),1);
    % out(:,i) = profileLikelihood(P(:,i),vals,stds,objfn_constants,profile_params,save_all_pars);
end
fprintf("FevalQueue finished.\n")
%%
for i = 1:n_abm_vecs

    [idx,temp] = fetchNext(FF);
    out(:,idx) = temp;

    if mod(i,100)==0
        temp = toc(t_start);
        fprintf("Finished %d after %s. ETR: %s\n",i,duration(0,0,temp),duration(0,0,temp/i * (n_abm_vecs-i)))
    end

end