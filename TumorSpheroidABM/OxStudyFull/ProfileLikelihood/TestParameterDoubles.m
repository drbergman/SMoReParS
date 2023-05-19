clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.reprint = true;


load("data/Profiles_SMFromData_clean.mat","out")
sm_par_display_names = ["\lambda","\alpha","K","d_{G1/S}","d_{G2/M}","EC50"];
sm_par_file_names = ["lambda","alpha","K","dG1S","dG2M","ec50"];

figure_layout = "individual";

colors = lines(2); % for when I need to use this to specifiy a color

switch figure_layout
    case "unified" % plot all triples on single figure
        f = figureOnRight("Name","ParDoublesAll","Units","pixels","Position",[0 0 1440 820]);
        tiledlayout("flow")
        np = size(out,1);
        for xi = 1:np
            for yi = (xi+1):np
                nexttile; hold on;
                % scatter(out{xi}(xi,:),out{xi}(yi,:),"filled","DisplayName",sm_par_display_names(xi))
                % scatter(out{yi}(xi,:),out{yi}(yi,:),"filled","DisplayName",sm_par_display_names(yi))
                plot(out{xi}(xi,:),out{xi}(yi,:),"Marker","o","MarkerFaceColor",colors(1,:),"MarkerSize",4,"DisplayName",sm_par_display_names(xi))
                plot(out{yi}(xi,:),out{yi}(yi,:),"Marker","o","MarkerFaceColor",colors(2,:),"MarkerSize",4,"DisplayName",sm_par_display_names(yi))
                xlabel(sm_par_display_names(xi))
                ylabel(sm_par_display_names(yi))
                L = legend("Location","best","FontSize",16);
                title(L,"Profiled Parameter")
            end
        end

    case "individual" % plot all triples on own figure
        save_fig_opts.subfolder = "ParameterDoubles";
        np = size(out,1);
        fi = 0;
        for xi = 1:np
            for yi = (xi+1):np
                fi = fi+1;
                f(fi) = figureOnRight("Name",sprintf("ParDouble_%s_%s",sm_par_file_names(xi),sm_par_file_names(yi)));
                hold on;
                % scatter(out{xi}(xi,:),out{xi}(yi,:),"filled","DisplayName",sm_par_display_names(xi))
                % scatter(out{yi}(xi,:),out{yi}(yi,:),"filled","DisplayName",sm_par_display_names(yi))
                plot(out{xi}(xi,:),out{xi}(yi,:),"Marker","o","MarkerFaceColor",colors(1,:),"MarkerSize",6,"DisplayName",sm_par_display_names(xi))
                plot(out{yi}(xi,:),out{yi}(yi,:),"Marker","o","MarkerFaceColor",colors(2,:),"MarkerSize",6,"DisplayName",sm_par_display_names(yi))
                xlabel(sm_par_display_names(xi))
                ylabel(sm_par_display_names(yi))
                L = legend("Location","best","FontSize",16);
                title(L,"Profiled Parameter")
                set(gca,"FontSize",16)
            end
        end
end

%% save figures
saveFigures(f,save_fig_opts);
