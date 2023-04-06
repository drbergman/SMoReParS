% plots histograms of the sensitivities of the parameters in both models

clearvars;

ABM = load("data/ABM_Sensitivity_noapop.mat");
ABM_SM = load("data/ABM_Sensitivity_Via_SM_noapop_UPDATED.mat");
figure;
bar(categorical(regexprep(ABM.par_names(ABM.order),"_"," "),regexprep(ABM.par_names(ABM.order),"_"," ")),[ABM.mu_star;ABM_SM.mu_star])
ylabel("\mu^*")
legend({"Sensitivity on ABM directly","Sensitivity on ABM via SM"},"location","northeast","FontSize",16)
set(gca,'FontSize',16)
savefig("figs/CompareABMSensitivityMethods")
print("figs/CompareABMSensitivityMethods","-dpng")
% load("data/ODE_Sensitivity_noapop.mat")
% figure;
% bar(categorical(regexprep(par_names(order),"_"," "),regexprep(par_names(order),"_"," ")),mu_star)