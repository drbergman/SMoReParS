clearvars;

save_opts.save_figs = true;
save_opts.reprint = true;
save_opts.file_types = ["fig","png"];
save_opts.fig_names = "SampleProfilesOfSMFromData_LMS";

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

opts = struct();
% opts.abm_vec_inds = 1:3; % set this to fix a subset of ABM parameters to test profiles
load("../ODEFitting/data/SMFitToData_LMS.mat","fixed_pars","model_type");
D = parameterDisplayNameDictionary(model_type);
sm_par_display_names = ["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"low_dose_apop";"delta_dose_apop";"rho0"];

profile_file = "data/Profiles_SMFromData_LMS_clean.mat";
nsamps = 1;

%% set up parameter names
[~,I] = setdiff(sm_par_display_names,fixed_pars);
sm_par_display_names = sm_par_display_names(sort(I));
for i = 1:numel(sm_par_display_names)
    sm_par_display_names(i) = D(sm_par_display_names(i));
end

%% test the profile
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);

%% save the figures
saveFigures(f,save_opts)

%% reset the path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../ODEFitting/")

