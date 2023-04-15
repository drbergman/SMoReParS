% plots histograms of the sensitivities of the parameters in both models

clearvars;

load("data/ABM_Sensitivity_noapop.mat")
figure;
bar(categorical(regexprep(par_names(order),"_"," "),regexprep(par_names(order),"_"," ")),mu_star)

load("data/ODE_Sensitivity_noapop.mat")
figure;
bar(categorical(regexprep(par_names(order),"_"," "),regexprep(par_names(order),"_"," ")),mu_star)