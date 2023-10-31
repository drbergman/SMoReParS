% this script identify admissible points in ABM parameter space from SMoRe 
% ParS.

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../../../SampleABMFns/")

% opts.acceptance_method = "single_best";
% opts.acceptance_method = "all_profiles";
% opts.acceptance_method = "all_profiles_resampled";
opts.acceptance_method = "all_profiles_resampled_enforce_ci_bounds";
% opts.par_combos = {1,2,3};
opts.nsamples = 10000;

files.optimal_parameters = "../ODEFitting/data/SMFitToABM_New.mat";

load(files.optimal_parameters,"cohort_name")

files.abm_data = sprintf("../../data/%s/summary_short.mat",cohort_name);
files.profile_from_abm = "../ProfileLikelihood/data/Profiles_SMFromABM_New_clean.mat";
files.profile_from_data = "../ProfileLikelihood/data/Profiles_SMFromData_New_clean.mat";

%% accept samples
accepted_parameters = acceptSampledABMParameters(files,opts);

%% save output
suffix = fileSuffix(opts);
file_name = "data/AcceptedParameters_New_" + suffix + ".mat";

save(file_name,"opts","cohort_name","files","accepted_parameters")

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../SampleABMFns/")

