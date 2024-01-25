clearvars;

addpath("../../ProfileLikelihoodFns/")

% model_type = "exponential";
% model_type = "logistic";
model_type = "von_bertalanffy";

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = "SampleProfilesOfSMFromABM_" + model_type;
save_fig_opts.reprint = true;

opts.abm_vec_inds = [51,57,35];
switch model_type
    case "exponential"
        sm_par_display_names = "\lambda";
    case "logistic"
        sm_par_display_names = ["r","K"];
    case "von_bertalanffy"
        sm_par_display_names = ["\alpha","\nu","\beta"];
end
profile_file = sprintf("data/ProfileLikelihoods_%s.mat",model_type);
nsamps = 3;
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);

saveFigures(f,save_fig_opts)


rmpath("../../ProfileLikelihoodFns/")

