clearvars

save_figs = true;
reprint = true;
file_types = ["fig";"png"];
resolution = "-r300";

resample_vb = true;

model_type = ["exponential","logistic","von_bertalanffy"];
sm_model_palette = smModelPalette();
f = figureOnRight("Name","RSSComparisonAcrossSMs");
ax = gca;
hold on
RSS = cell(3,1);
for i = 1:3
    if model_type(i)~="von_bertalanffy" || ~resample_vb
        temp = load(sprintf("data/OptimalParameters_%s.mat",model_type(i)));
    else
        temp = load(sprintf("data/OptimalParameters_%s_resampled.mat",model_type(i)));
    end
    RSS{i} = temp.fstar;
    if isfield(temp,"resample_t")
        RSS_bar = RSS{i} ./ length(temp.resample_t);
    else
        RSS_bar = RSS{i} ./ 300;
    end
    h = histogram(log10(RSS{i}(:) / length(temp.resample_t)),"FaceColor",sm_model_palette(model_type(i)));
    if i==3
        h.morebins
        h.morebins
        h.morebins
    end
end
xlabel("log_{10}(RSS)")
ylabel("Frequency")
xticks(-4:2:2)
display_names = ["Exp","Log","vB"];
set(ax,"FontSize",8)
% legend(display_names)

%%
f.Units = "inches";
f.Position(3:4) = [2,1.5];

%%
saveFigures(f,save_figs=save_figs, reprint=reprint, file_types=file_types, resolution=resolution)