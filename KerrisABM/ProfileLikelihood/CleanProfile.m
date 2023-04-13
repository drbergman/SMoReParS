clearvars;

addpath("../../ProfileLikelihoodFns/")
load("data/ProfileLikelihoods.mat","out")
threshold = chi2inv(0.95,3);

out = cleanProfiles(out,threshold);

save("data/ProfileLikelihoods.mat","out")

rmpath("../../ProfileLikelihoodFns/")

