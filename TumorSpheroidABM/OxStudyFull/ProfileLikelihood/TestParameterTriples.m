clearvars;

save_fig_opts.save_figs = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.reprint = true;

file_base_name = "Data";

load("../ODEFitting/data/SMFitTo" + file_base_name,"fixed_pars")
load("data/Profiles_SMFrom" + file_base_name + "_clean.mat","profiles")
sm_par_display_names = ["\lambda";"\alpha";"K";"\alpha_R";"\alpha_P";"k_\alpha";"a";"\delta_0";"k_\delta";"b";"\rho_0"];
for i = 1:numel(fixed_pars)
    sm_par_display_names(regexprep(sm_par_display_names,'\','')==fixed_pars(i)) = [];
end
sm_par_file_names = regexprep(sm_par_display_names,'\','');

figure_layout = "individual";

colors = lines(3); % for when I need to use this to specifiy a color

switch figure_layout
    case "unified" % plot all triples on single figure
        f = figureOnRight("Name","ParTriplesAll","Units","pixels","Position",[0 0 1440 820]);
        tiledlayout("flow")

        np = size(profiles,1);
        for xi = 1:np
            for yi = (xi+1):np
                for zi = (yi+1):np
                    nexttile; hold on;
                    % scatter3(profiles{xi}(xi,:),profiles{xi}(yi,:),profiles{xi}(zi,:),"filled","DisplayName",sm_par_display_names(xi))
                    % scatter3(profiles{yi}(xi,:),profiles{yi}(yi,:),profiles{yi}(zi,:),"filled","DisplayName",sm_par_display_names(yi))
                    % scatter3(profiles{zi}(xi,:),profiles{zi}(yi,:),profiles{zi}(zi,:),"filled","DisplayName",sm_par_display_names(zi))
                    plot3(profiles{xi}(xi,:),profiles{xi}(yi,:),profiles{xi}(zi,:),"Marker","o","MarkerFaceColor",colors(1,:),"MarkerSize",2,"DisplayName",sm_par_display_names(xi))
                    plot3(profiles{yi}(xi,:),profiles{yi}(yi,:),profiles{yi}(zi,:),"Marker","o","MarkerFaceColor",colors(2,:),"MarkerSize",2,"DisplayName",sm_par_display_names(yi))
                    plot3(profiles{zi}(xi,:),profiles{zi}(yi,:),profiles{zi}(zi,:),"Marker","o","MarkerFaceColor",colors(3,:),"MarkerSize",2,"DisplayName",sm_par_display_names(zi))
                    xlabel(sm_par_display_names(xi))
                    ylabel(sm_par_display_names(yi))
                    zlabel(sm_par_display_names(zi))
                    view(3)
                    L = legend("Location","best","FontSize",12);
                    title(L,"Profiled Parameter")
                end
            end
        end

    case "individual" % plot all triples on own figure
        save_fig_opts.subfolder = "ParameterTriples";
        np = size(profiles,1);
        fi = 0;
        for xi = 1:np
            for yi = (xi+1):np
                for zi = (yi+1):np
                    fi = fi+1;
                    f(fi) = figureOnRight("Name",sprintf("ParTriple_%s_%s_%s",sm_par_file_names(xi),sm_par_file_names(yi),sm_par_file_names(zi)));
                    hold on;
                    % scatter3(profiles{xi}(xi,:),profiles{xi}(yi,:),profiles{xi}(zi,:),"filled","DisplayName",sm_par_display_names(xi))
                    % scatter3(profiles{yi}(xi,:),profiles{yi}(yi,:),profiles{yi}(zi,:),"filled","DisplayName",sm_par_display_names(yi))
                    % scatter3(profiles{zi}(xi,:),profiles{zi}(yi,:),profiles{zi}(zi,:),"filled","DisplayName",sm_par_display_names(zi))
                    plot3(profiles{xi}(xi,:),profiles{xi}(yi,:),profiles{xi}(zi,:),"Marker","o","MarkerFaceColor",colors(1,:),"MarkerSize",2,"DisplayName",sm_par_display_names(xi))
                    plot3(profiles{yi}(xi,:),profiles{yi}(yi,:),profiles{yi}(zi,:),"Marker","o","MarkerFaceColor",colors(2,:),"MarkerSize",2,"DisplayName",sm_par_display_names(yi))
                    plot3(profiles{zi}(xi,:),profiles{zi}(yi,:),profiles{zi}(zi,:),"Marker","o","MarkerFaceColor",colors(3,:),"MarkerSize",2,"DisplayName",sm_par_display_names(zi))
                    xlabel(sm_par_display_names(xi))
                    ylabel(sm_par_display_names(yi))
                    zlabel(sm_par_display_names(zi))
                    view(3)
                    L = legend("Location","best","FontSize",16);
                    title(L,"Profiled Parameter")
                    set(gca,"FontSize",16)
                end
            end
        end
end

%% save figures
saveFigures(f,save_fig_opts);
