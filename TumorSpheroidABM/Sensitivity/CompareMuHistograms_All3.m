% plots histograms of the sensitivities of the parameters in both models

clearvars;

ABM = load("ABM_Sensitivity_noapop.mat");
ABM_SM = load("ABM_Sensitivity_Via_SM_noapop.mat");
figure;
bar(categorical(regexprep(ABM.par_names(ABM.order),"_"," "),regexprep(ABM.par_names(ABM.order),"_"," ")),[ABM.mu_star;ABM_SM.mu_star])

load("ODE_Sensitivity_noapop.mat")
figure;
bar(categorical(regexprep(par_names(order),"_"," "),regexprep(par_names(order),"_"," ")),mu_star)