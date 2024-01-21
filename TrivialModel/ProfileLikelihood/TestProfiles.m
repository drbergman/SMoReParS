clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../ProfileLikelihoodFns/")

profile_file = "data/Profiles_clean.mat";

nsamps = 12;
sm_par_display_names = "s";
opts = struct();

[f,ax,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);

rmpath("../../ProfileLikelihoodFns/")
