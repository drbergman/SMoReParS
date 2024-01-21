clearvars;

save_figure = false;

addpath("../../ProfileLikelihoodFns/")

model_type = "exponential";

switch model_type
    case "exponential"
        sm_par_display_names = "\lambda";
    case "logistic"
        sm_par_display_names = ["r","K"];
    case "von_bertalanffy"
        sm_par_display_names = ["\alpha","\nu","\beta"];
end
profile_file = sprintf("data/ProfileLikelihoods_%s.mat",model_type);
n_per_fig = 10;
f = testAllProfileSMFromABM(profile_file,n_per_fig,sm_par_display_names);

if save_figure
    for i = 1:numel(f)
        savefig(f,sprintf("figures/fig/AllProfilesOfSMFromABM_%d",i))
        print(f,sprintf("figures/png/SampleProfilesOfSMFromABM_%d",i),"-dpng")
    end
end

rmpath("../../ProfileLikelihoodFns/")

