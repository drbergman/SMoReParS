clearvars

addpath("~/Documents/MATLAB/myfunctions/")

file_name = "SMFitToData_New";
experimental_data = "data/ExperimentalData_New.mat";

save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = file_name;
save_fig_opts.resolution = '-r1200';

fit_color = lines(2);
fit_color = sqrt(prod(fit_color,1));
%% prepare stuff
load("data/" + file_name,"P","fn_opts")
load(experimental_data,"t","D");

%%
f=figure;
hold on;
tfull = linspace(0,3,100);
yy = computeTimeSeries(P,tfull,[],fn_opts,zeros(6,1));
data_curve = plot(t,D.A,"black","Marker","o","MarkerFaceColor","black","MarkerSize",3,"DisplayName","Data","LineStyle","--");
data_patch = patch([t;flip(t)],[D.A-D.S;flip(D.A+D.S)],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
fit_curve = plot(tfull,sum(yy,2),"Color",fit_color,"LineStyle","-","LineWidth",0.5,"DisplayName","Fit");
xlabel("Time (d)")
ylabel("Cell Count")

% legend([fit_curve;data_curve;data_patch],"Location","northwest","FontSize",8)

f.Units = "inches";
f.Position(3) = 2;
f.Position(4) = 1;
set(gca,"FontSize",8)

%% adjust margin
margin = struct("left",.18,"right",.03,"top",[],"bottom",.28);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(gca,margin,spacing);

%% save figure
saveFigures(f,save_fig_opts)