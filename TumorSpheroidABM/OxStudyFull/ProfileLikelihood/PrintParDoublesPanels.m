clearvars;

% this script will print paper-ready versions of the parameter doubles

save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.subfolder = "ParameterDoubles_SMFitToData_LMS_bounded";
save_fig_opts.file_types = "png"; % leave the original alone
save_fig_opts.resolution = "-r1200";

par_pairs = ["lambda","alpha";"lambda","K";"alpha","K"; % control parameters
             "alphaR","kalpha";"alphaR","a";"alphaR","kdelta";"kdelta","delta"];

all_figs = dir("figures/fig/ParameterDoubles_SMFitToData_LMS_bounded/*.fig");
all_fig_names = {all_figs.name};
n = size(par_pairs,1);
f = gobjects(n,1);

for i = 1:n
    contains_first = contains(all_fig_names,"_" + par_pairs(i,1) + "_") | contains(all_fig_names,"_" + par_pairs(i,1) + ".");
    contains_second = contains(all_fig_names,"_" + par_pairs(i,2) + "_") | contains(all_fig_names,"_" + par_pairs(i,2) + ".");
    file_log = contains_first & contains_second;
    file_name = all_fig_names{file_log};

    f(i) = open("figures/fig/ParameterDoubles_SMFitToData_LMS_bounded/" + file_name);

    delete(f(i).Children(1))
    ax = f(i).Children;
    set(ax,"FontSize",8);

    f(i).Units = "inches";
    f(i).Position(3) = 1.5;
    f(i).Position(4) = 1;

    if ax.YLim(2) == 10000
        ax.YAxis.Exponent = 3;
    end

    set(ax.Children,"MarkerSize",4,"Marker",".")

    %% adjust margin
    margin = struct("left",.21,"right",.07,"top",.15,"bottom",.34);
    spacing = struct("horizontal",0.1,"vertical",0.19);
    uniformAxisSpacing(ax,margin,spacing);

end

%%
saveFigures(f,save_fig_opts)