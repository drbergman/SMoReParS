clearvars;

addpath("../../../ProfileLikelihoodFns/")
% load("data/Profiles_SMFromABM.mat","out") % profiles from ABM
load("data/Profiles_SMFromData.mat","out") % profiles from data
% load("data/Profiles_SMFromABM_extended.mat","out") % profiles from data
threshold = chi2inv(0.95,6);

out = cleanProfiles(out,threshold);

save("data/Profiles_SMFromData.mat","out")
% save("data/Profiles_SMFromData_clean.mat","out")

rmpath("../../../ProfileLikelihoodFns/")
