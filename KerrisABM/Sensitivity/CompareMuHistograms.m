% plots histograms of the sensitivities of the parameters in both models

clearvars;

save_fig_opts.file_types = ["fig","png"];

ABM = load("data/GlobalSensitivityMOATDirect.mat");
ABM_SM = load("data/GlobalSensitivityMOATIndirect_logistic.mat");

C = categorical(ABM.ordered_par_names,ABM.ordered_par_names);
line_width = 1;

%% get indirect order to match direct order
[~,order_sm] = sort(ABM_SM.display_par_names);
[~,order_abm] = sort(ABM.ordered_par_names);
order_abm_inv(order_abm) = 1:length(order_abm);
temp = ABM_SM.mu_star(:)';
temp = temp(order_sm); % put it alphabetical order
temp = temp(order_abm_inv); % put it in same order as direct method
y = [ABM.mu_star(:),temp'];

%% raw mu
f=figureOnRight("Name","CompareABMSensitivityMethods");
hold on;
b=bar(y,"LineWidth",line_width);
xticks(1:size(y,1))
xticklabels(C)
[ngroups,nbars] = size(y);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
std_error = [ABM.sigma(:),ABM_SM.sigma(:)] ./ [sqrt(ABM.npoints),sqrt(ABM_SM.npoints)];
er = errorbar(x',y,zeros(size(y)),std_error,"CapSize",20);    
set(er,"Color",[0 0 0],"LineStyle","none","LineWidth",line_width);                            
ylabel("\mu^*")
legend({"Sensitivity on ABM directly","Sensitivity on ABM via SM"},"location","northeast","FontSize",16)
set(gca,'FontSize',16)

%% normalized mu
y_sum = sum(y,1); 
y_normalized = y./y_sum;
f(2) = figureOnRight("Name","CompareABMSensitivityMethods_NormalizedMu");
hold on;
b=bar(y_normalized,"LineWidth",line_width);
xticks(1:size(y_normalized,1))
xticklabels(C)
[ngroups,nbars] = size(y_normalized);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
std_error = ([ABM.sigma(:),ABM_SM.sigma(:)] ./ [sqrt(ABM.npoints),sqrt(ABM_SM.npoints)])./y_sum;
er = errorbar(x',y_normalized,zeros(size(y_normalized)),std_error,"CapSize",20);    
set(er,"Color",[0 0 0],"LineStyle","none","LineWidth",line_width); 
ylabel("\mu^* / <\mu^*>")
legend({"Sensitivity on ABM directly","Sensitivity on ABM via SM"},"location","northeast","FontSize",16)
set(gca,'FontSize',16)

%% stacked histogram
f(3) = figureOnRight("Name","CompareABMSensitivityMethods_NormalizedMu_Stacked");
s = ["Direct","Indirect"];
B = bar(categorical(s,s),y_normalized,"stacked");
ylabel("\mu^* / <\mu^*>")
legend(flip(B),flip(C),"location","bestoutside","FontSize",16)
set(gca,"FontSize",16)

%% save figs
saveFigures(f,save_fig_opts)
