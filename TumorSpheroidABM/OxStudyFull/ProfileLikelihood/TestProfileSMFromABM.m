clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

profile_file = "data/Profiles_SMFromABM_LMS_clean.mat";

save_opts.save_figs = true;
save_opts.reprint = true;
save_opts.file_types = ["fig","png"];
save_opts.fig_names = "SampleProfilesOfSMFromABM_LMS";

load("../ODEFitting/data/SMFitToData_LMS","fixed_pars","model_type")

opts = struct();
% opts.abm_vec_inds = 1:3; % set this to fix a subset of ABM parameters to test profiles
sm_par_display_names = ["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"low_dose_apop";"delta_dose_apop";"rho0"];

nsamps = 5;

%% set up parameter names
D = parameterDisplayNameDictionary(model_type);
[~,I] = setdiff(sm_par_display_names,fixed_pars);
sm_par_display_names = sm_par_display_names(sort(I));
for i = 1:numel(sm_par_display_names)
    sm_par_display_names(i) = D(sm_par_display_names(i));
end

%% make plots
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);

%% save figures
saveFigures(f,save_opts)

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../ODEFitting/")

