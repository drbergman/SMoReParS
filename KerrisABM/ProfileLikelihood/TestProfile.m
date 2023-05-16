clearvars;

save_figure = false;
is_cleaned = false;

addpath("../../ProfileLikelihoodFns/")

model_type = "von_bertalanffy";

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = sprintf("SampleProfilesOfSMFromABM_%s",model_type);

switch model_type
    case "logistic"
        sm_par_display_names = ["r","K"];
    case "von_bertalanffy"
        sm_par_display_names = ["\alpha","\nu","\beta"];
end
profile_file = sprintf("data/ProfileLikelihoods_%s.mat",model_type);
nsamps = 5;
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names);

saveFigures(f,save_fig_opts)


rmpath("../../ProfileLikelihoodFns/")

