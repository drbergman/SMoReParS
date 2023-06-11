% this script identify admissible points in ABM parameter space from SMoRe 
% ParS.

clearvars;

addpath("../../../ProfileLikelihoodFns/")
addpath("../../../SampleABMFns/")


% opts.admission_method = "single_best";
opts.admission_method = "all_profiles_resampled";
opts.nsamples = 500; % converges to 2040 admitted by setting this to 750 (but 700 is not big enough)

files.par_file = "../ODEFitting/data/SMFitToABM_LMS.mat";

load(files.par_file,"cohort_name")

files.abm_data_file = sprintf("../../data/%s/summary.mat",cohort_name);
files.profile_from_abm_file = "../ProfileLikelihood/data/Profiles_SMFromABM_LMS_clean.mat";
files.profile_from_data_file = "../ProfileLikelihood/data/Profiles_SMFromData_LMS_clean.mat";

%% admit samples
admitted_parameters = admitSampledABMParameters(files,opts);

%% save output
% save("data/AdmittedParameters_LMS_" + opts.admission_method + ".mat","opts","cohort_name","files","admitted_parameters")

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../SampleABMFns/")


