% plots histograms of the sensitivities of the parameters in both models

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../ODEFitting/")
addpath("..")

save_figs = true;
reprint = true;
file_types = ["fig","png"];
resolution = "-r300";

direct_method_color = 0.1*[1 1 1];
sm_model_palette = smModelPalette();

% model_type = "exponential";
% model_type = "logistic";
model_type = ["exponential","logistic","von_bertalanffy"];
include_vB_in_normalized = false; % no longer show vB in normalized figs since it's so bad
legend_entries = ["Direct","Exp","Log","vB"];
% suffix = "";
% suffix = "_large";
suffix = "_very_large";

% endpoint = "final_size";
% endpoint = "AUC";
endpoint = "time_to_half";

if endpoint == "time_to_half"
    model_type(model_type=="von_bertalanffy") = [];
    legend_entries(legend_entries=="vB") = [];
end

ABM = load(sprintf("data/GlobalSensitivityMOATDirect_%s.mat",endpoint));
for i = 1:length(model_type)
    if model_type(i)=="von_bertalanffy"
        ABM_SM(i) = load(sprintf("data/GlobalSensitivityMOATIndirect_%s%s_%s_resampled_clean.mat",model_type(i),suffix,endpoint),"display_par_names","mu_star","sigma","npoints");
    else
        ABM_SM(i) = load(sprintf("data/GlobalSensitivityMOATIndirect_%s%s_%s.mat",model_type(i),suffix,endpoint),"display_par_names","mu_star","sigma","npoints");
    end
end

C = categorical(ABM.ordered_par_names,ABM.ordered_par_names);
line_width = 0.5;

%% get indirect order to match direct order
y = zeros(size(ABM.mu_star,1),1+length(model_type));
y(:,1) = ABM.mu_star(:);
for i = 1:length(model_type)
    [~,order_sm] = sort(ABM_SM(i).display_par_names);
    [~,order_abm] = sort(ABM.ordered_par_names);
    order_abm_inv = zeros(1,length(order_abm));
    order_abm_inv(order_abm) = 1:length(order_abm);
    temp = ABM_SM(i).mu_star(:)';
    temp = temp(order_sm); % put it alphabetical order
    temp = temp(order_abm_inv); % put it in same order as direct method
    y(:,i+1) = temp;
end

%% raw mu
f=figureOnRight("Name",sprintf("CompareMOAT%s_%s",suffix,endpoint));
ax = gca;
hold on;
% y_log = log10(y);
b=bar(y,"LineWidth",line_width);
% b=bar(y_log,"LineWidth",line_width);
b(1).FaceColor = direct_method_color;
for i = 1:length(model_type)
    b(i+1).FaceColor = sm_model_palette(model_type(i));
end
% set(b,"EdgeColor","none")
xticks(1:size(y,1))
xticklabels(C)
[ngroups,nbars] = size(y);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
std_error = zeros(size(y));
std_error(:,1) = ABM.sigma ./ sqrt(ABM.npoints);
for i = 1:length(model_type)
    std_error(:,i+1) = ABM_SM(i).sigma ./ sqrt(ABM_SM(i).npoints);
end
% std_error_log = log10(std_error);
% std_error = [ABM.sigma(:),ABM_SM.sigma(:)] ./ [sqrt(ABM.npoints),sqrt(ABM_SM.npoints)];
er = errorbar(x',y,zeros(size(y)),std_error,"CapSize",10*0.257142857142857);    
% er = errorbar(x',y_log,zeros(size(y_log)),std_error_log,"CapSize",10);    
% er = errorbar(x',y_log,zeros(size(y_log)),log10(y+std_error)-y_log,"CapSize",10);    
set(er,"Color",[0 0 0],"LineStyle","none","LineWidth",line_width);                            
% ylabel("\mu^*")
% legend(legend_entries,"location","northeast","FontSize",16)
set(gca,'FontSize',6)
% yT = 0:10:50;
% yticks(yT)
% for i = 1:length(yT)
%     yTL(i) = sprintf("10^{%d}",yT(i));
% end
% yticklabels(yTL)
set(ax,"YScale","log")
drawnow
if endpoint == "final_size"
    yticks(10.^[150,200])
elseif endpoint == "AUC"
    yticks(10.^[8,10])
    exp_10 = (10.^[8,10]) .* log10(exp(1));
end

f.Units = "inches";
f.Position(3:4) = [2,1.5];

saveFigures(f,save_figs=save_figs, reprint=reprint, resolution=resolution, file_types=file_types)
set(ax,"YScale","linear")
if endpoint == "final_size"
    b(4).YData(:) = 1e7; % prevents weird behavior in vB bars on this zoomed in scale
    ylim([0 2e5])
    yticks([0 1e5])
    yticklabels(["0","10^5"]);
elseif endpoint == "AUC"
    ylim([0 1e7])
    yticks([0 5e6])
    yticklabels(["0","5x10^6"]);
end


f.Name = f.Name + "_zoom";
saveFigures(f,save_figs=save_figs, reprint=reprint, resolution=resolution, file_types=file_types)


%% normalized mu
y_sum = sum(y,1); 
y_normalized = y./y_sum;
if ~include_vB_in_normalized && any(model_type=="von_bertalanffy")
    y_normalized(:,end) = [];
end
% f = figureOnRight("Name",sprintf("CompareABMSensitivityMethods_NormalizedMu_%s%s",model_type,suffix));
% hold on;
% b=bar(y_normalized,"LineWidth",line_width);
% xticks(1:size(y_normalized,1))
% xticklabels(C)
% [ngroups,nbars] = size(y_normalized);
% x = nan(nbars, ngroups);
% for i = 1:nbars
%     x(i,:) = b(i).XEndPoints;
% end
% std_error = ([ABM.sigma(:),ABM_SM.sigma(:)] ./ [sqrt(ABM.npoints),sqrt(ABM_SM.npoints)])./y_sum;
% er = errorbar(x',y_normalized,zeros(size(y_normalized)),std_error,"CapSize",20);    
% set(er,"Color",[0 0 0],"LineStyle","none","LineWidth",line_width); 
% ylabel("\mu^* / <\mu^*>")
% legend({"Sensitivity on ABM directly","Sensitivity on ABM via SM"},"location","northeast","FontSize",16)
% set(gca,'FontSize',16)

%% stacked histogram
f = figureOnRight("Name",sprintf("CompareABMSensitivityMethods_NormalizedMu_Stacked_All%s_%s_legend",suffix,endpoint));
stacked_hist_legend = legend_entries;
if ~include_vB_in_normalized
    stacked_hist_legend(stacked_hist_legend=="vB") = [];
end
B = bar(categorical(stacked_hist_legend,stacked_hist_legend),y_normalized',"stacked");
% ylabel("\mu^* / <\mu^*>")
legend(flip(B),flip(C),"location","bestoutside","FontSize",8)
abm_par_palette = abmParameterPalette();
for i = 1:numel(B)
    B(i).FaceColor = abm_par_palette(B(i).DisplayName);
end
set(gca,"FontSize",16)
ylim([0 1])
yticks([0 1])
% legend("off")

%%
ax = gca;
ax.Legend.FontSize = 6;

f.Units = "inches";
f.Position(3:4) = [2,1.5];
set(gca,"FontSize",6)

%% save figs
saveFigures(f,save_figs=save_figs, reprint=reprint, resolution=resolution, file_types=file_types)

%% save without legend
f.Name = sprintf("CompareABMSensitivityMethods_NormalizedMu_Stacked_All%s_%s",suffix,endpoint);
legend off
saveFigures(f,save_figs=save_figs, reprint=reprint, resolution=resolution, file_types=file_types)

%% reset path
rmpath("../ODEFitting/")
rmpath("..")

