% A script to run the indirect global sensitivity.

clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../ODEFitting/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../../../SensitivityFns/")

use_profiles = true;
sort_output = false;

Nr = 15; % number of resamples per factor in ABM space
nsamps = 200; % number of points to sample in LHS for ODE pars
omega_max = 8;
M = 4;
Ns = 65;
% Ns = 249;

cohort_name = "cohort_230124175743017";
files.profiles = "../ProfileLikelihood/data/Profiles_SMFromABM_New_clean.mat";

% PL = load("../ProfileLikelihood/data/Profiles_SMFromABM_New_clean.mat","profiles");
% load(sprintf("../../data/%s/summary.mat",cohort_name),"vals","cohort_size","par_names")
% par_names = ["carrying_capacity";"occmax_2d";"move_rate_microns";"g1_to_s";"s_to_g2";"g2_to_m";"m_to_g1"];

D = makeABMParameterDistributionsDictionary();
T = dictionary("occmax_2d",@(x) min(7,floor(x)));

sm_functional = @(p) sum(computeTimeSeries(p, [], [], false, 3));

%% run MOAT
[studied_function,par_names] = setupSampleFromSMFunction(files,sm_functional,T=T,D=D,nsamps=nsamps,use_profiles=use_profiles);
n_abm_pars = length(par_names);
[S1,ST] = efast(studied_function,n_abm_pars,Nr,omega_max,M,Ns);

%% save result
save("data/GlobalSensitivityeFASTIndirect2.mat","S1","ST",...
    "par_names","Nr","nsamps","Ns","M","omega_max")

%% clean path
rmpath("../ODEFitting/")
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../SensitivityFns/")


