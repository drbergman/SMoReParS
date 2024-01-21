clearvars

addpath("../../ProfileLikelihoodFns/")


overwrite_profile = false;
profile_to_clean = "data/Profiles";

%% load and clean P.profiles
P = load(profile_to_clean); % P.profiles from ABM
threshold = chi2inv(0.95,size(P.profiles,1));

P.profiles = cleanProfiles(P.profiles,threshold);

if overwrite_profile
    save(profile_to_clean,"-struct","P") %#ok<*UNRCH>
else
    save(profile_to_clean + "_clean","-struct","P")
end

%% reset path
rmpath("../../ProfileLikelihoodFns/")
