clearvars;

save_figure = false;
is_cleaned = true;

addpath("../../ProfileLikelihoodFns/")


sm_par_display_names = ["\alpha","\theta","\beta"];
profile_file = "data/ProfileLikelihoods.mat";
nsamps = 10;
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,is_cleaned);

if save_figure
    savefig(f,"figures/fig/SampleProfilesOfSMFromABM")
    print(f,"figures/png/SampleProfilesOfSMFromABM","-dpng")
end

rmpath("../../ProfileLikelihoodFns/")

