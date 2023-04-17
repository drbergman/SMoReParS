clearvars;

save_figure = false;

addpath("../../../ProfileLikelihoodFns/")

sm_par_display_names = ["\lambda","\alpha","K"];
profile_file = "data/ProfileLikelihoods.mat";
n_per_fig = 10;
f = testAllProfileSMFromABM(profile_file,n_per_fig,sm_par_display_names);

if save_figure
    for i = 1:numel(f)
        savefig(f,sprintf("figures/fig/AllProfilesOfSMFromABM_%d",i))
        print(f,sprintf("figures/png/SampleProfilesOfSMFromABM_%d",i),"-dpng")
    end
end

rmpath("../../ProfileLikelihoodFns/")

