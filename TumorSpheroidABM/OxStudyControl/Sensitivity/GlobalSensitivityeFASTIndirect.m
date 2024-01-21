% A script to run the indirect global sensitivity.

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../ODEFitting/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../../../SensitivityFns/")


Nr = 15; % number of resamples per factor in ABM space
nsamps = 200; % number of points to sample in LHS for ODE pars
omega_max = 8;
M = 4;
Ns = 65;
% Ns = 249;

cohort_name = "cohort_230124175743017";

PL = load("../ProfileLikelihood/data/Profiles_SMFromABM_New_clean.mat","profiles");
load(sprintf("../../data/%s/summary.mat",cohort_name),"vals","cohort_size","par_names")

n_abm_pars = length(par_names);
D = makeABMParameterDistributionsDictionary(par_names);
T = dictionary("occmax_2d",@(x) min(7,floor(x)));

% end_fn = @(v) v(end);
sm_functional = @(x) sum(computeTimeSeries(x, [], [], false, 3));
%% create bounding surfaces
n_sm_pars = size(PL.profiles,1);
PL.profiles = reshape(PL.profiles,n_sm_pars,[]);
n_abm_vecs = size(PL.profiles,2);

BS = zeros(n_sm_pars,n_abm_vecs,2);
threshold = chi2inv(0.95,n_sm_pars);
for i = 1:n_abm_vecs
    for j = 1:n_sm_pars
        [BS(j,i,1),BS(j,i,2)] = getProfileBounds(PL.profiles{j,i}([j,end],:),threshold);
    end
end
BS = reshape(BS,[n_sm_pars,cohort_size,2]);

%% run MOAT
studied_function = @(x) sampleFromSM(x,BS,vals,sm_functional,D=D,T=T,nsamps=nsamps,par_names=par_names);
% studied_function = @(x) moatSample_ODE(x,par_names,D);
[S1,ST,ST_desc_order] = efast(studied_function,n_abm_pars,Nr,omega_max,M,Ns);
[~,S1_desc_order] = sort(S1,"descend");
par_names_desc_S1_order = par_names(S1_desc_order);
par_names_desc_ST_order = par_names(ST_desc_order);

%% save result
save("data/GlobalSensitivityeFASTIndirect.mat","S1","ST","S1_desc_order",...
    "ST_desc_order","par_names","par_names_desc_S1_order","par_names_desc_ST_order","Nr","nsamps","Ns","M","omega_max")

%% clean path
rmpath("../ODEFitting/")
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../SensitivityFns/")


