clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.reprint = true;


load("data/Profiles_SMFromData_clean.mat","out")
sm_par_display_names = ["\lambda","\alpha","K","d_{G1/S}","d_{G2/M}","EC50"];
sm_par_file_names = ["lambda","alpha","K","dG1S","dG2M","ec50"];

figure_layout = "individual";

colors = lines(3); % for when I need to use this to specifiy a color

switch figure_layout
    case "unified" % plot all triples on single figure
        f = figureOnRight("Name","ParTriplesAll","Units","pixels","Position",[0 0 1440 820]);
        tiledlayout("flow")

        np = size(out,1);
        for xi = 1:np
            for yi = (xi+1):np
                for zi = (yi+1):np
                    nexttile; hold on;
                    % scatter3(out{xi}(xi,:),out{xi}(yi,:),out{xi}(zi,:),"filled","DisplayName",sm_par_display_names(xi))
                    % scatter3(out{yi}(xi,:),out{yi}(yi,:),out{yi}(zi,:),"filled","DisplayName",sm_par_display_names(yi))
                    % scatter3(out{zi}(xi,:),out{zi}(yi,:),out{zi}(zi,:),"filled","DisplayName",sm_par_display_names(zi))
                    plot3(out{xi}(xi,:),out{xi}(yi,:),out{xi}(zi,:),"Marker","o","MarkerFaceColor",colors(1,:),"MarkerSize",2,"DisplayName",sm_par_display_names(xi))
                    plot3(out{yi}(xi,:),out{yi}(yi,:),out{yi}(zi,:),"Marker","o","MarkerFaceColor",colors(2,:),"MarkerSize",2,"DisplayName",sm_par_display_names(yi))
                    plot3(out{zi}(xi,:),out{zi}(yi,:),out{zi}(zi,:),"Marker","o","MarkerFaceColor",colors(3,:),"MarkerSize",2,"DisplayName",sm_par_display_names(zi))
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
        np = size(out,1);
        fi = 0;
        for xi = 1:np
            for yi = (xi+1):np
                for zi = (yi+1):np
                    fi = fi+1;
                    f(fi) = figureOnRight("Name",sprintf("ParTriple_%s_%s_%s",sm_par_file_names(xi),sm_par_file_names(yi),sm_par_file_names(zi)));
                    hold on;
                    % scatter3(out{xi}(xi,:),out{xi}(yi,:),out{xi}(zi,:),"filled","DisplayName",sm_par_display_names(xi))
                    % scatter3(out{yi}(xi,:),out{yi}(yi,:),out{yi}(zi,:),"filled","DisplayName",sm_par_display_names(yi))
                    % scatter3(out{zi}(xi,:),out{zi}(yi,:),out{zi}(zi,:),"filled","DisplayName",sm_par_display_names(zi))
                    plot3(out{xi}(xi,:),out{xi}(yi,:),out{xi}(zi,:),"Marker","o","MarkerFaceColor",colors(1,:),"MarkerSize",2,"DisplayName",sm_par_display_names(xi))
                    plot3(out{yi}(xi,:),out{yi}(yi,:),out{yi}(zi,:),"Marker","o","MarkerFaceColor",colors(2,:),"MarkerSize",2,"DisplayName",sm_par_display_names(yi))
                    plot3(out{zi}(xi,:),out{zi}(yi,:),out{zi}(zi,:),"Marker","o","MarkerFaceColor",colors(3,:),"MarkerSize",2,"DisplayName",sm_par_display_names(zi))
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
