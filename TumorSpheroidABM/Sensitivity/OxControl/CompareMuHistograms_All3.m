% plots histograms of the sensitivities of the parameters in both models

clearvars;

ABM = load("data/ABM_Sensitivity_noapop.mat");
ABM_SM = load("data/ABM_Sensitivity_Via_SM_noapop_UPDATED.mat");

C = regexprep(ABM.par_names(ABM.order),"_"," ");
C = categorical(C,C);
y = [ABM.mu_star;ABM_SM.mu_star];

figure;
bar(C,y)
ylabel("\mu^*")
legend({"Sensitivity on ABM directly","Sensitivity on ABM via SM"},"location","northeast","FontSize",16)
set(gca,'FontSize',16)
savefig("figures/fig/CompareABMSensitivityMethods")
print("figures/png/CompareABMSensitivityMethods","-dpng")

y = y./sum(y,2);
figure;
bar(C,y)
ylabel("\mu^* / <\mu^*>")
legend({"Sensitivity on ABM directly","Sensitivity on ABM via SM"},"location","northeast","FontSize",16)
set(gca,'FontSize',16)
savefig("figures/fig/CompareABMSensitivityMethods_NormalizedMu")
print("figures/png/CompareABMSensitivityMethods_NormalizedMu","-dpng")

%%
figure;
s = ["Direct","Indirect"];
B = bar(categorical(s,s),y,"stacked");
ylabel("\mu^* / <\mu^*>")
legend(flip(B),flip(C),"location","bestoutside","FontSize",16)
set(gca,"FontSize",16)
savefig("figures/fig/CompareABMSensitivityMethods_NormalizedMu_Stacked")
print("figures/png/CompareABMSensitivityMethods_NormalizedMu_Stacked","-dpng")
