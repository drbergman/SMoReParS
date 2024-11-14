clearvars

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../../SurrogateModelFns/")
addpath("..")

clear rawError

force_serial = true;
n_starts = 1;

cohort_name = "cohort_1";

sm.custom_solve_sm_fn = @(sm,p,tt,C,D,condition_on_previous,resample_t) surrogateModel(p);
p = 0;
lb = -Inf;
ub = Inf;

files.data = sprintf("../data/%s/summary.mat",cohort_name);
optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
weights = 1;

P = optimizeSMParsFromABM(files,sm,p,lb,ub,optim_opts,weights,...
    force_serial=force_serial,n_starts=n_starts);

save("data/OptimalParameters.mat","P","files","sm","p","lb","ub","optim_opts","weights","n_starts")

rmpath("../../SurrogateModelFns/")
rmpath("..")