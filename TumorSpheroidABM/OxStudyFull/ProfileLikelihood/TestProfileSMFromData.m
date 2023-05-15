clearvars;

save_opts.save_figs = true;
save_opts.reprint = true;
save_opts.file_types = ["fig","png"];
save_opts.fig_names = "SampleProfilesOfSMFromData";

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ProfileLikelihoodFns/")

opts = struct();
% opts.abm_vec_inds = 1:3; % set this to fix a subset of ABM parameters to test profiles
sm_par_display_names = ["\lambda","\alpha","K","d_{G1/S}","d_{G2/M}","EC50"];
profile_file = "data/Profiles_SMFromData_clean.mat";
nsamps = 1;
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);

saveFigures(f,save_opts)

rmpath("../../../ProfileLikelihoodFns/")

