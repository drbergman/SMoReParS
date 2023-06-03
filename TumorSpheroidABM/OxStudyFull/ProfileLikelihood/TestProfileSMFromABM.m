clearvars;

profile_file = "data/Profiles_SMFromABM_clean.mat";

save_opts.save_figs = true;
save_opts.reprint = true;
save_opts.file_types = ["fig","png"];
save_opts.fig_names = "SampleProfilesOfSMFromABM";

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ProfileLikelihoodFns/")

opts = struct();
% opts.abm_vec_inds = 1:3; % set this to fix a subset of ABM parameters to test profiles
sm_par_display_names = ["\alpha_R";"\alpha_P";"k_\alpha";"\delta_0";"k_\delta";"\rho_0"];
nsamps = 5;
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);


saveFigures(f,save_opts)

rmpath("../../../ProfileLikelihoodFns/")

