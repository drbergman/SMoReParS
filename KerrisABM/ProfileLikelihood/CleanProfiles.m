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

P = load(file_name);
threshold = chi2inv(0.95,npars);

P.profiles = cleanProfiles(P.profiles,threshold);

save(file_name,"-struct","P")

rmpath("../../ProfileLikelihoodFns/")

