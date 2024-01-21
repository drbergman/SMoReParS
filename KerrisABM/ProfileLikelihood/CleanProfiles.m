clearvars;

addpath("../../ProfileLikelihoodFns/")

model_type = "exponential";
file_name = sprintf("data/ProfileLikelihoods_%s",model_type);
switch model_type
    case "exponential"
        npars = 1;
    case "logistic"
        npars = 2;
    case "von_bertalanffy"
        npars = 3;
end

load(file_name,"profiles")
threshold = chi2inv(0.95,npars);

profiles = cleanProfiles(profiles,threshold);

save(file_name,"profiles")

rmpath("../../ProfileLikelihoodFns/")

