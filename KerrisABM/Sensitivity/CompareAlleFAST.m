% plots histograms of the sensitivities of the parameters in both models

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../ODEFitting/")
addpath("..")

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.resolution = "-r1200";

direct_method_color = 0.1*[1 1 1];
sm_model_palette = smModelPalette();

% model_type = "exponential";
% model_type = "logistic";
model_type = ["exponential","logistic","von_bertalanffy"];
legend_entries = ["Direct","Exponential","Logistic","Von Bertalanffy"];
% suffix = "";
% suffix = "_large";
suffix = "_big_sample";

% endpoint = "final_size";
% endpoint = "AUC";
endpoint = "time_to_half";

if endpoint == "time_to_half"
    model_type(model_type=="von_bertalanffy") = [];
    legend_entries(legend_entries=="Von Bertalanffy") = [];
end

ABM = load(sprintf("data/GlobalSensitivityeFASTDirect_%s.mat",endpoint));
for i = 1:length(model_type)
    ABM_SM(i) = load(sprintf("data/GlobalSensitivityeFASTIndirect_%s%s_%s.mat",model_type(i),suffix,endpoint),"display_par_names","S1","ST");
end

line_width = 1;

%% get indirect order to match direct order
y_S1 = zeros(size(ABM.S1,1),1+length(model_type));
[y_S1(:,1),order_abm] = sort(ABM.S1,"descend");
C = categorical(ABM.display_par_names(order_abm),ABM.display_par_names(order_abm));
for i = 1:length(model_type)
    y_S1(:,i+1) = ABM_SM(i).S1(order_abm);
end

%% raw mu
f=figureOnRight("Name",sprintf("CompareeFAST_S1_%s_%s",suffix,endpoint));
hold on;
b=bar(y_S1,"LineWidth",line_width);
b(1).FaceColor = direct_method_color;
for i = 1:length(model_type)
    b(i+1).FaceColor = sm_model_palette(model_type(i));
end
xticks(1:size(y_S1,1))
xticklabels(C)
set(gca,'FontSize',16)
ylim([0 1])
yticks([0 1])
saveFigures(f,save_fig_opts)

%% get indirect order to match direct order
y_ST = zeros(size(ABM.ST,1),1+length(model_type));
[y_ST(:,1),order_abm] = sort(ABM.ST,"descend");
C = categorical(ABM.display_par_names(order_abm),ABM.display_par_names(order_abm));
for i = 1:length(model_type)
    y_ST(:,i+1) = ABM_SM(i).ST(order_abm);
end

%% raw mu
f=figureOnRight("Name",sprintf("CompareeFAST_ST_%s_%s",suffix,endpoint));
hold on;
b=bar(y_ST,"LineWidth",line_width);
b(1).FaceColor = direct_method_color;
for i = 1:length(model_type)
    b(i+1).FaceColor = sm_model_palette(model_type(i));
end
xticks(1:size(y_ST,1))
xticklabels(C)
[ngroups,nbars] = size(y_ST);
set(gca,'FontSize',16)
ylim([0 1])
yticks([0 1])
saveFigures(f,save_fig_opts)


%% normalized mu
y_S1_sum = sum(y_S1,1);
y_S1_normalized = y_S1./y_S1_sum;

%% stacked histogram
f = figureOnRight("Name",sprintf("CompareeFAST_S1_Methods_NormalizedMu_Stacked_All%s_%s",suffix,endpoint));
B = bar(categorical(legend_entries,legend_entries),y_S1_normalized',"stacked");
% ylabel("\mu^* / <\mu^*>")
legend(flip(B),flip(C),"location","bestoutside","FontSize",16)
abm_par_palette = abmParameterPalette();
for i = 1:numel(B)
    B(i).FaceColor = abm_par_palette(B(i).DisplayName);
end
set(gca,"FontSize",16)
ylim([0 1])
yticks([0 1])

%% save figs
saveFigures(f,save_fig_opts)

%% normalized mu
y_ST_sum = sum(y_ST,1);
y_ST_normalized = y_ST./y_ST_sum;

%% stacked histogram
f = figureOnRight("Name",sprintf("CompareeFAST_ST_Methods_NormalizedMu_Stacked_All%s_%s",suffix,endpoint));
B = bar(categorical(legend_entries,legend_entries),y_ST_normalized',"stacked");
% ylabel("\mu^* / <\mu^*>")
legend(flip(B),flip(C),"location","bestoutside","FontSize",16)
abm_par_palette = abmParameterPalette();
for i = 1:numel(B)
    B(i).FaceColor = abm_par_palette(B(i).DisplayName);
end
set(gca,"FontSize",16)
ylim([0 1])
yticks([0 1])

%% save figs
saveFigures(f,save_fig_opts)

%% reset path
rmpath("../ODEFitting/")
rmpath("..")

