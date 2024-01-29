clearvars

addpath("~/Documents/MATLAB/myfunctions/")

save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig";"png"];
save_fig_opts.resolution = "-r1200";

model_type = ["exponential","logistic","von_bertalanffy"];

identifiability_index_palette = identifiabilityIndexPalette();
color_order = [identifiability_index_palette("0");identifiability_index_palette("1");identifiability_index_palette("2")];
f = gobjects(6,1);
k = 0;
for i = 1:length(model_type)
    load(sprintf("data/IdentifiabilityIndex_%s",model_type(i)),"indices")
    for j = 1:size(indices,1)
        k = k+1;
        f(k) = figureOnRight("Name",sprintf("IdentifiabilityIndex_%s_%d",model_type(i),j));
        C = categorical(indices(j,:),0:2);
        d = donutchart(C,InnerRadius=0.2);
        d.Labels(:) = "";
        % d.Labels(d.CategoryCounts==0) = "";
        colororder(color_order)
    end
end

%%
scale_factor = 4;
set(f,"Units","inches")
for i = 1:numel(f)
    f(i).Position(3:4) = 0.5 * scale_factor;
    f(i).Children.FontSize = 6 * scale_factor;
end

%%
saveFigures(f,save_fig_opts)
