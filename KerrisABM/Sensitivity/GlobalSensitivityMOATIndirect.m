% A script to run the indirect global sensitivity.

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../ODEFitting/")
addpath("../../ProfileLikelihoodFns/")
addpath("../../SensitivityFns/")

npoints = 1000; % number of points to sample in LHS for ABM pars
nsamps = 100; % number of points to sample in LHS for ODE pars

suffix = dictionary([15,25,1000],["","_large","_very_large"]);

% model_type = "exponential";
% model_type = "logistic";
model_type = "von_bertalanffy";


sm_functional = @(p) computeSMEndpoint(p,[],model_type);
switch model_type
    case "exponential"
        sum_fn = @mean;
    case "logistic"
        sum_fn = @mean; % use this for logistic growth
    case "von_bertalanffy"
        sum_fn = @summarizeSMLHS; % use this for the VB model when the parameters make it likely that simulations blow up
end

files.profiles = sprintf("../ProfileLikelihood/data/ProfileLikelihoods_%s.mat",model_type);

% PL = load(sprintf("../ProfileLikelihood/data/ProfileLikelihoods_%s.mat",model_type),"profiles");
load("../PostAnalysis/data/summary.mat","vals","cohort_size","display_par_names")

n_abm_pars = length(display_par_names);
D = makeMOATDistributions(display_par_names);
T = makeParameterTransformations(display_par_names);

%% create bounding surfaces
% n_sm_pars = size(PL.profiles,1);
% PL.profiles = reshape(PL.profiles,n_sm_pars,[]);
% n_abm_vecs = size(PL.profiles,2);
% 
% BS = zeros(n_sm_pars,n_abm_vecs,2);
% threshold = chi2inv(0.95,n_sm_pars);
% for i = 1:n_abm_vecs
%     for j = 1:n_sm_pars
%         [BS(j,i,1),BS(j,i,2)] = getProfileBounds(PL.profiles{j,i}([j,end],:),threshold);
%     end
% end
% BS = reshape(BS,[n_sm_pars,cohort_size,2]);

%% run MOAT
% studied_function = @(x) sampleFromSM(x,BS,vals,sm_functional,D=D,T=T,nsamps=nsamps,par_names=display_par_names,sum_fn=sum_fn);
% studied_function = @(x) sampleFromSM(x,files,sm_functional,D=D,T=T,nsamps=nsamps,sum_fn=sum_fn);
studied_function = setupSampleFromSMFunction(files,sm_functional,D=D,T=T,nsamps=nsamps,sum_fn=sum_fn,par_names = display_par_names, warnings=false);
[mu_star,sigma,order] = morris_simple(studied_function,n_abm_pars,npoints);
display_par_names = display_par_names(order);

%% save result
save(sprintf("data/GlobalSensitivityMOATIndirect_%s%s.mat",model_type,suffix(npoints)),"mu_star","sigma","display_par_names","npoints","nsamps")

%% clean path
rmpath("../ODEFitting/")
rmpath("../../ProfileLikelihoodFns/")
rmpath("../../SensitivityFns/")


