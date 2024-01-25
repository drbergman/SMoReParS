clearvars

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig";"png"];
save_fig_opts.resolution = "-r1200";

model_type = ["exponential","logistic","von_bertalanffy"];
sm_model_palette = smModelPalette();
f = figureOnRight("Name","RSSComparisonAcrossSMs");
hold on

for i = 1:3
    temp = load(sprintf("data/OptimalParameters_%s.mat",model_type(i)));
    fstar = temp.fstar;
    if isfield(temp,"resample_t")
        fstar = fstar ./ length(temp.resample_t);
    else
        fstar = fstar ./ 300;
    end
    histogram(log10(fstar(:)),"FaceColor",sm_model_palette(model_type(i)))
end

display_names = ["Exponential","Logistic","Von Bertalanffy"];
legend(display_names)

saveFigures(f,save_fig_opts)