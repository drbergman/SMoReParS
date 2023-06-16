clearvars;

addpath("../../../ProfileLikelihoodFns/")

overwrite_profile = false;
profile_to_clean = "data/Profiles_SMFromABM_New";

%% load and clean profiles
load(profile_to_clean,"profiles") % profiles from ABM
threshold = chi2inv(0.95,size(profiles,1));

profiles = cleanProfiles(profiles,threshold);

if overwrite_profile
    save(profile_to_clean,"profiles") %#ok<*UNRCH>
else
    save(profile_to_clean + "_clean","profiles")
end

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
