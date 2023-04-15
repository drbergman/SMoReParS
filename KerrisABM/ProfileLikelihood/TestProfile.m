clearvars;

save_figure = true;
is_cleaned = true;

addpath("../../ProfileLikelihoodFns/")


sm_par_display_names = ["\alpha","\nu","\beta"];
profile_file = "data/ProfileLikelihoods.mat";
nsamps = 5;
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names);

if save_figure
    savefig(f,"figures/fig/SampleProfilesOfSMFromABM")
    print(f,"figures/png/SampleProfilesOfSMFromABM","-dpng")
end

rmpath("../../ProfileLikelihoodFns/")

