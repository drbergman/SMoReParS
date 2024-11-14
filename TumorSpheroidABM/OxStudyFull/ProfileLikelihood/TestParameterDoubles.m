clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../ODEFitting/")

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.reprint = true;
% save_fig_opts.resolution = '-r1200';

show_legend = true;

file_base_name = "Data_LMS_bounded";

load("../ODEFitting/data/SMFitTo" + file_base_name,"fixed_pars","model_type")
load("data/Profiles_SMFrom" + file_base_name + "_clean_converted.mat","profiles")
sm_par_file_names = ["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"delta";"kdelta";"rho0"];

figure_layout = "individual"; % unified or individual

colors = lines(2); % for when I need to use this to specifiy a color


%% set up parameter names
D = parameterDisplayNameDictionary(model_type);
[~,I] = setdiff(sm_par_file_names,fixed_pars);
sm_par_file_names = sm_par_file_names(sort(I));
sm_par_display_names = sm_par_file_names;
for i = 1:numel(sm_par_display_names)
    sm_par_display_names(i) = D(sm_par_display_names(i));
end

%% plots
switch figure_layout
    case "unified" % plot all doubles on single figure
        f = figureOnRight("Name","ParDoublesAll_SMFitTo" + file_base_name,"Units","pixels","Position",[0 0 1440 820]);
        tiledlayout("flow")
        np = size(profiles,1);
        for xi = 1:np
            for yi = (xi+1):np
                nexttile; hold on;
                % scatter(profiles{xi}(xi,:),profiles{xi}(yi,:),"filled","DisplayName",sm_par_display_names(xi))
                % scatter(profiles{yi}(xi,:),profiles{yi}(yi,:),"filled","DisplayName",sm_par_display_names(yi))
                plot(profiles{xi}(xi,:),profiles{xi}(yi,:),"Marker","o","MarkerFaceColor",colors(1,:),"MarkerSize",4,"DisplayName",sm_par_display_names(xi))
                plot(profiles{yi}(xi,:),profiles{yi}(yi,:),"Marker","o","MarkerFaceColor",colors(2,:),"MarkerSize",4,"DisplayName",sm_par_display_names(yi))
                xlabel(sm_par_display_names(xi),"FontSize",16)
                ylabel(sm_par_display_names(yi),"FontSize",16)
                if show_legend
                    L = legend("Location","best","FontSize",16);
                    title(L,"Profiled Parameter")
                end
            end
        end
        if ~show_legend
            L = legend(["x","y"],"Location","best","FontSize",16);
            title(L,"Profiled Dimension")
            fprintf("You probably want to move the legend off the final plot. Pausing for you to do that...\n")
            pause
        else
            fprintf("You probably want to move legends off the plots.\n")
        end
    case "individual" % plot all doubles on own figure
        save_fig_opts.subfolder = "ParameterDoubles_SMFitTo" + file_base_name;
        np = size(profiles,1);
        fi = 0;
        for xi = 1:np
            for yi = (xi+1):np
                fi = fi+1;
                f(fi) = figureOnRight("Name",sprintf("ParDouble_%s_%s",sm_par_file_names(xi),sm_par_file_names(yi)));
                hold on;
                % scatter(profiles{xi}(xi,:),profiles{xi}(yi,:),"filled","DisplayName",sm_par_display_names(xi))
                % scatter(profiles{yi}(xi,:),profiles{yi}(yi,:),"filled","DisplayName",sm_par_display_names(yi))
                plot(profiles{xi}(xi,:),profiles{xi}(yi,:),"Marker","o","MarkerFaceColor",colors(1,:),"MarkerSize",6,"DisplayName",sm_par_display_names(xi))
                plot(profiles{yi}(xi,:),profiles{yi}(yi,:),"Marker","o","MarkerFaceColor",colors(2,:),"MarkerSize",6,"DisplayName",sm_par_display_names(yi))
                xlabel(sm_par_display_names(xi))
                ylabel(sm_par_display_names(yi))
                if show_legend
                    L = legend("Location","best","FontSize",16);
                    title(L,"Profiled Parameter")
                end
                set(gca,"FontSize",16)
            end
        end
end

%% save figures
saveFigures(f,save_fig_opts);

%% reset path
rmpath("../ODEFitting/")
