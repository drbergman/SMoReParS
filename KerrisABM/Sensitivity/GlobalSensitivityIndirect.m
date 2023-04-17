% A script to run the indirect global sensitivity.

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../ODEFitting/")
addpath("../../ProfileLikelihoodFns/")
addpath("../../SensitivityFns/")

npoints = 15; % number of points to sample in LHS for ABM pars
nsamps = 100; % number of points to sample in LHS for ODE pars

fn = @computeSMEndpoint;
sum_fn = @summarizeSMLHS;

PL = load("../ProfileLikelihood/data/ProfileLikelihoods.mat","out");
load("../PostAnalysis/summary.mat","vals","cohort_size","display_par_names")

n_abm_pars = length(display_par_names);
D = makeMOATDistributions(display_par_names);
T = makeParameterTransformations(display_par_names);

%% create bounding surfaces
n_sm_pars = size(PL.out,1);
PL.out = reshape(PL.out,n_sm_pars,[]);
n_abm_vecs = size(PL.out,2);

BS = zeros(n_sm_pars,n_abm_vecs,2);
threshold = chi2inv(0.95,3);
for i = 1:n_abm_vecs
    for j = 1:n_sm_pars
        [BS(j,i,1),BS(j,i,2)] = getProfileBounds(PL.out{j,i}([j,end],:),threshold);
    end
end
BS = reshape(BS,[n_sm_pars,cohort_size,2]);

%% run MOAT
studied_function = @(x) moatSampleFromSM(x,display_par_names,BS,T,D,vals,nsamps,fn,{[]},[],sum_fn);
[mu_star,sigma,order] = morris_simple(studied_function,n_abm_pars,npoints);
display_par_names = display_par_names(order);

%% save result
save("data/GlobalSensitivityIndirect.mat","mu_star","sigma","display_par_names","npoints")

%% clean path
rmpath("../ODEFitting/")
rmpath("../../ProfileLikelihoodFns/")
rmpath("../../SensitivityFns/")


