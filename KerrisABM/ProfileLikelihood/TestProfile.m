clearvars;

addpath("../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

model_type = "exponential";
% model_type = "logistic";
% model_type = "von_bertalanffy";

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = "SampleProfilesOfSMFromABM_" + model_type;
save_fig_opts.reprint = true;

sm_color_palette = smModelPalette();

opts.abm_vec_inds = [51,57,35];
opts.LineWidth = 0.5;
opts.place_par_names = "none";
opts.LineColor = sm_color_palette(model_type);
% opts.show_y_label = false;
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
[f,ax] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);

nc = size(ax,2);
horizontal_margin = 0.2;
horizontal_space = 1.25*horizontal_margin;
panel_width = 0.5;
scale_factor = 1;
f.Units = "inches";
% f.Position(3) = (2*horizontal_margin + horizontal_space*(nc-1) + panel_width * nc) * scale_factor;
% fig_width = nc * panel_width / (1- 2*horizontal_margin + (nc-1)*horizontal_space) * scale_factor;
fig_width = (2*horizontal_margin + (nc-1)*horizontal_space + nc * panel_width) * scale_factor;
f.Position(3) = fig_width;
f.Position(4) = 2 * scale_factor;
set(f.Children,"FontSize" ,6 * scale_factor);

rel_horizontal_margin = horizontal_margin/fig_width;
rel_horizontal_space = horizontal_space / fig_width;
%% margins for sample fit figure
margin = struct("left",rel_horizontal_margin,"right",rel_horizontal_margin,"top",.05,"bottom",0.15);
spacing = struct("horizontal",rel_horizontal_space,"vertical",0.15);
uniformAxisSpacing(ax,margin,spacing);


%% save figure
saveFigures(f,save_fig_opts)


rmpath("../../ProfileLikelihoodFns/")
rmpath("../ODEFitting/")

