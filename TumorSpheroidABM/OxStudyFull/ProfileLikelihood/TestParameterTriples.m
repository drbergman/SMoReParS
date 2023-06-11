clearvars;

addpath("../ODEFitting/")

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.reprint = true;

show_legend = true;

file_base_name = "Data_LMS";

load("../ODEFitting/data/SMFitTo" + file_base_name,"fixed_pars","model_type")
load("data/Profiles_SMFrom" + file_base_name + "_clean.mat","profiles")
sm_par_file_names = ["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"low_dose_apop";"delta_dose_apop";"rho0"];

figure_layout = "individual"; % unified or individual

colors = lines(3); % for when I need to use this to specifiy a color

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
    case "unified" % plot all triples on single figure
        f = figureOnRight("Name","ParTriplesAll","Units","pixels","Position",[0 0 1440 820]);
        tiles = tiledlayout("flow");

        np = size(profiles,1);
        for xi = 1:np
            for yi = (xi+1):np
                for zi = (yi+1):np
                    nexttile; hold on;plot3(profiles{xi}(xi,:),profiles{xi}(yi,:),profiles{xi}(zi,:),"Marker","o","MarkerFaceColor",colors(1,:),"MarkerSize",2,"DisplayName",sm_par_display_names(xi))
                    plot3(profiles{yi}(xi,:),profiles{yi}(yi,:),profiles{yi}(zi,:),"Marker","o","MarkerFaceColor",colors(2,:),"MarkerSize",2,"DisplayName",sm_par_display_names(yi))
                    plot3(profiles{zi}(xi,:),profiles{zi}(yi,:),profiles{zi}(zi,:),"Marker","o","MarkerFaceColor",colors(3,:),"MarkerSize",2,"DisplayName",sm_par_display_names(zi))
                    xlabel(sm_par_display_names(xi))
                    ylabel(sm_par_display_names(yi))
                    zlabel(sm_par_display_names(zi))
                    view(3)
                    if show_legend
                        L = legend("Location","best","FontSize",12);
                        title(L,"Profiled Parameter")
                    end
                end
            end
        end
        if ~show_legend
            for i = 1:numel(tiles.Children)
                xlabel(tiles.Children(i),sprintf("x=%s",tiles.Children(i).XLabel.String))
                ylabel(tiles.Children(i),sprintf("y=%s",tiles.Children(i).YLabel.String))
                zlabel(tiles.Children(i),sprintf("z=%s",tiles.Children(i).ZLabel.String))
            end
            L = legend(["x","y","z"],"Location","best","FontSize",16);
            title(L,"Profiled Dimension")
            fprintf("You probably want to move the legend off the final plot. Pausing for you to do that...\n")
            pause
        end

    case "individual" % plot all triples on own figure
        save_fig_opts.subfolder = "ParameterTriples";
        np = size(profiles,1);
        f = gobjects(1,nchoosek(np,3));
        fi = 0;
        for xi = 1:np
            for yi = (xi+1):np
                for zi = (yi+1):np
                    fi = fi+1;
                    f(fi) = figureOnRight("Name",sprintf("ParTriple_%s_%s_%s",sm_par_file_names(xi),sm_par_file_names(yi),sm_par_file_names(zi)));
                    hold on;plot3(profiles{xi}(xi,:),profiles{xi}(yi,:),profiles{xi}(zi,:),"Marker","o","MarkerFaceColor",colors(1,:),"MarkerSize",2,"DisplayName",sm_par_display_names(xi))
                    plot3(profiles{yi}(xi,:),profiles{yi}(yi,:),profiles{yi}(zi,:),"Marker","o","MarkerFaceColor",colors(2,:),"MarkerSize",2,"DisplayName",sm_par_display_names(yi))
                    plot3(profiles{zi}(xi,:),profiles{zi}(yi,:),profiles{zi}(zi,:),"Marker","o","MarkerFaceColor",colors(3,:),"MarkerSize",2,"DisplayName",sm_par_display_names(zi))
                    view(3)
                    xlabel(sm_par_display_names(xi))
                    ylabel(sm_par_display_names(yi))
                    zlabel(sm_par_display_names(zi))
                    if show_legend
                        L = legend("Location","best","FontSize",16);
                        title(L,"Profiled Parameter")
                    end
                    set(gca,"FontSize",16)
                end
            end
        end
end

%% save figures
saveFigures(f,save_fig_opts);

%% reset path
rmpath("../ODEFitting/")