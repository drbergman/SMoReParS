clearvars;

save_opts.save_figs = false;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ProfileLikelihoodFns/")

opts = struct();
% opts.abm_vec_inds = 1:3; % set this to fix a subset of ABM parameters to test profiles
sm_par_display_names = ["\lambda","\alpha","K","d_{G1/S}","d_{G2/M}","EC50"];
% profile_file = "data/Profiles_SMFromABM_old.mat";
profile_file = "data/Profiles_SMFromABM.mat";
nsamps = 10;
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);


save_opts.file_types = ["fig","png"];
save_opts.fig_names = "SampleProfilesOfSMFromABM";
saveFigures(f,save_opts)

rmpath("../../../ProfileLikelihoodFns/")

