clearvars;

% this script will test out how to do sensitivity on the CM using the SM

%% Program to run

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../ProfileLikelihoodFns/")
addpath("../../SensitivityFns/")
addpath("..")

options.limit_factor = 0.5; % how to set the limit for separating low and high impact factors
options.initialization_factor = 0.8; % how to determine when it has been sufficiently initialized
nsim_max = 1e4;
nsamps = 100; % number of points to sample in LHS for ODE pars
use_profiles = false;

par_names = ["a";"b"];
[D,I] = makeCMParameterDistributionsDictionary(par_names,distribution="Normal");

%% create bounding hypersurfaces
% cohort_name = "cohort_1";
files.profiles = "../ProfileLikelihood/data/Profiles_clean.mat";
% PL = load("../ProfileLikelihood/data/Profiles_clean.mat");
% C = load(sprintf("../data/%s/summary.mat",cohort_name),"cohort_size","vals");
% vals = C.vals;
% 
% npars_ode = size(PL.profiles,1);
% PL.profiles = reshape(PL.profiles,npars_ode,[]);
% npoints = size(PL.profiles,2);


% Number of factors of uncertainty of the function studied :
nfac=numel(par_names); 

sm_functional = @complexModel;
% studied_function = setupSampleFromSMFunction(files,sm_functional,par_names=par_names,D=D,nsamps=nsamps,use_profiles=use_profiles);
[studied_function,par_names] = setupSampleFromSMFunction(files,sm_functional,D=D,nsamps=nsamps,use_profiles=use_profiles);
[mu_star,sigma,order] = morris_simple(studied_function,nfac,15);

rmpath("..")
rmpath("../../ProfileLikelihoodFns/")
rmpath("../../SensitivityFns/")

