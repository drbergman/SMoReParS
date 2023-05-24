% this script identify admissible points in ABM parameter space from SMoRe 
% ParS.

clearvars;

addpath("../../../ProfileLikelihoodFns/")
addpath("../../../SampleABMFns/")


% opts.admission_method = "single_best";
opts.admission_method = "all_profiles_resampled";
opts.nsamples = 750; % converges to 2040 admitted by setting this to 750 (but 700 is not big enough)

cohort_name = "cohort_2303301105";

files.abm_data_file = sprintf("../../data/%s/summary_new.mat",cohort_name);
files.profile_from_abm_file = "../ProfileLikelihood/data/Profiles_SMFromABM_clean.mat";
files.profile_from_data_file = "../ProfileLikelihood/data/Profiles_SMFromData_clean.mat";

%% admit samples
admitted_parameters = admitSampledABMParameters(files,opts);

%% save output
save("data/AdmittedParameters_" + opts.admission_method + ".mat","opts","cohort_name","files","admitted_parameters")

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../SampleABMFns/")


