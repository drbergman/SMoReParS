clearvars

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../../../ProfileLikelihoodFns/")

save_figs = true;
reprint = true;
file_types = ["fig";"png"];
resolution = "-r300";

load("data/IdentifiabilityIndex","indices")
par_names = ["lambda","alpha","K"];

identifiability_index_palette = identifiabilityIndexPalette();
color_order = [identifiability_index_palette("0");identifiability_index_palette("1");identifiability_index_palette("2")];
n_sm_pars = size(indices,1);
f = gobjects(n_sm_pars,1);
k = 0;
for j = 1:n_sm_pars
    k = k+1;
    f(k) = figureOnRight("Name",sprintf("IdentifiabilityIndex_%s",par_names(j)));
    C = categorical(indices(j,:),0:2);
    d = donutchart(C,InnerRadius=0.6);
    d.Labels(:) = "";
    % d.Labels(d.CategoryCounts==0) = "";
    colororder(color_order)
end

%%
scale_factor = 4;
set(f,"Units","inches")
for i = 1:numel(f)
    f(i).Position(3:4) = 0.5 * scale_factor;
    f(i).Children.FontSize = 6 * scale_factor;
end

%%
saveFigures(f,save_figs=save_figs,reprint=reprint,file_types=file_types,resolution=resolution)


%% reset path
rmpath("../../../ProfileLikelihoodFns/")
