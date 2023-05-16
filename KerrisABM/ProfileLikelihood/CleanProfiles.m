clearvars;

addpath("../../ProfileLikelihoodFns/")

model_type = "logistic";
file_name = sprintf("data/ProfileLikelihoods_%s",model_type);
switch model_type
    case "logistic"
        npars = 2;
    case "von_bertalanffy"
        npars = 3;
end

load(file_name,"out")
threshold = chi2inv(0.95,npars);

out = cleanProfiles(out,threshold);

save(file_name,"out")

rmpath("../../ProfileLikelihoodFns/")

