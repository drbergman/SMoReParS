% this script accepts points in ABM parameter space from SMoRe 
% ParS.

clearvars;

addpath("../../../ProfileLikelihoodFns/")
addpath("../../../SampleABMFns/")

name_suffix = "LMS_2";

% opts.acceptance_method = "single_best";
opts.acceptance_method = "specified_parameter_combinations_resampled";
% opts.acceptance_method = "all_profiles_resampled";
opts.nsamples = 10000; % converges to 2040 accepted by setting this to 750 (but 700 is not big enough)
opts.par_combos = num2cell(1:9);
files.optimal_parameters = "../ODEFitting/data/SMFitToABM_" + name_suffix + ".mat";

load(files.optimal_parameters,"cohort_name")

files.abm_data = sprintf("../../data/%s/summary.mat",cohort_name);
files.profile_from_abm = "../ProfileLikelihood/data/Profiles_SMFromABM_" + name_suffix + "_clean.mat";
files.profile_from_data = "../ProfileLikelihood/data/Profiles_SMFromData_LMS_bounded_clean.mat";

%% accept samples
accepted_parameters = acceptSampledABMParameters(files,opts);

%% save output
suffix = fileSuffix(opts);
file_name = "data/AcceptedParameters_" + name_suffix + "_" + suffix;

save(file_name,"opts","cohort_name","files","accepted_parameters")

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../SampleABMFns/")


