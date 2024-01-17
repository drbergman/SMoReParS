clearvars;

addpath("../../../SurrogateModelFns/")
addpath("~/Documents/MATLAB/myfunctions/")

file_name = "SMFitToData_LMS_bounded";
save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = file_name;
save_fig_opts.resolution = '-r1200';

load("data/" + file_name + ".mat","P");
load("data/" + file_name + ".mat","fn","fn_opts","sm")
if ~exists("sm","var")
    sm.fn = fn;
    sm.opts = fn_opts;
end
load("data/ExperimentalData.mat","t","D","C");

%% color
fit_color = lines(2);
fit_color = sqrt(prod(fit_color,1));
opts.LineColor = fit_color;

%%
f=figureOnRight;
tfull = linspace(0,3,100);
ax = gobjects(2,3);
for i = 1:3
    for j = 1:2
        ax(j,i) = subplot(2,3,r2c(2,3,[j,i])); hold on;
    end
end
for i = 1:3
    sim_data = fn(P,tfull,C{i},fn_opts);
    plot(ax(1,i),tfull,sim_data(:,1),"-","LineWidth",0.5,"DisplayName","Fit","Color",fit_color);
    plot(ax(2,i),tfull,sim_data(:,2),"-","LineWidth",0.5,"DisplayName","Fit","Color",fit_color);
    patch(ax(1,i),[t;flip(t)],[D(i).A(:,1)-D(i).S(:,1);flip(D(i).A(:,1)+D(i).S(:,1))],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
    plot(ax(1,i),t,D(i).A(:,1),"black","Marker","*","MarkerSize",3,"MarkerFaceColor","black","DisplayName","Data","LineWidth",0.5,"LineStyle","--");
    patch(ax(2,i),[t;flip(t)],[D(i).A(:,2)-D(i).S(:,2);flip(D(i).A(:,2)+D(i).S(:,2))],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
    plot(ax(2,i),t,D(i).A(:,2),"black","Marker","*","MarkerSize",3,"MarkerFaceColor","black","DisplayName","Data","LineWidth",0.5,"LineStyle","--");
    if i==1
        title(ax(1,i),"Control")
    else
        title(ax(1,i),sprintf("C = %3.2f\\muM",C{i}))
    end
end
xlabel(ax(2,:),"Time (d)")
ylabel(ax(1,1),"Cell Count","FontWeight","bold")
ylabel(ax(2,1),"G2/M Fraction","FontWeight","bold")
normalizeYLims(ax(2,:))

% legend([fit_curve;data_curve;data_patch],"Location","northwest","FontSize",22)

set(ax,"FontSize",8)

%%
f.Units = "inches";
f.Position(3) = 2.5;
f.Position(4) = 2;

%% adjust margin
margin = struct("left",.15,"right",.05,"top",.08,"bottom",.15);
spacing = struct("horizontal",0.1,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%% remove control g2/m fraction
delete(ax(2,1).Children)
ax(2,1).Color = f.Color;
ax(2,1).XColor = "none";
ax(2,1).YColor = "none";
ax(2,1).YAxis.Label.Color = "black";

%%
saveFigures(f,save_fig_opts)

%% remove paths
rmpath("../../../SurrogateModelFns/")
