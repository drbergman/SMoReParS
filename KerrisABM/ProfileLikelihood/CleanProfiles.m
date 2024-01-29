clearvars;

addpath("../../ProfileLikelihoodFns/")

model_type = "von_bertalanffy";
files.profiles = sprintf("data/ProfileLikelihoods_%s",model_type);
boundary_tolerance = [[1e-4;1e-4;1e-4],[1e-1;1;1e-1]];

overwrite_profile = false;
switch model_type
    case "exponential"
        npars = 1;
    case "logistic"
        npars = 2;
    case "von_bertalanffy"
        npars = 3;
end

% P = load(file_name);
% threshold = chi2inv(0.95,npars);

profiles = cleanProfiles(files, boundary_tolerance=boundary_tolerance);

if overwrite_profile
    file_name = files.profiles;
else
    file_name = files.profiles + "_clean";
    copyfile(files.profiles + ".mat",file_name + ".mat")
end

save(file_name,"profiles","boundary_tolerance","-append") %#ok<*UNRCH>

rmpath("../../ProfileLikelihoodFns/")

