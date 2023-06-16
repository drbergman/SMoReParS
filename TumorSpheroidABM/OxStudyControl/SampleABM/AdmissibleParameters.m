% this script identify admissible points in ABM parameter space from SMoRe 
% ParS.

clearvars;

addpath("../../../ProfileLikelihoodFns/")
addpath("../../../SampleABMFns/")

% opts.admission_method = "single_best";
% opts.admission_method = "all_profiles";
% opts.admission_method = "all_profiles_resampled";
opts.admission_method = "specified_parameter_combinations_resampled";
opts.par_combos = {[1,2]};
opts.nsamples = 10000; % converges to 2040 admitted by setting this to 750 (but 700 is not big enough)

files.optimal_parameters = "../ODEFitting/data/SMFitToABM_New.mat";

load(files.optimal_parameters,"cohort_name")

files.abm_data = sprintf("../../data/%s/summary_short.mat",cohort_name);
files.profile_from_abm = "../ProfileLikelihood/data/Profiles_SMFromABM_New_clean.mat";
files.profile_from_data = "../ProfileLikelihood/data/Profiles_SMFromData_New_clean.mat";

%% admit samples
admitted_parameters = admitSampledABMParameters(files,opts);

%% save output
file_name = "data/AdmittedParameters_New_" + opts.admission_method + ".mat";
save(file_name,"opts","cohort_name","files","admitted_parameters")
if isfield(opts,"admission_method") && opts.admission_method == "specified_parameter_combinations"
    par_combos = opts.par_combos;
    save(file_name,"par_combos","-append")
end

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../SampleABMFns/")


