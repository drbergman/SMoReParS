clearvars;

load("ABM_Sensitivity_noapop.mat")
figure;
bar(categorical(regexprep(par_names(order),"_"," "),regexprep(par_names(order),"_"," ")),mu_star)

load("ODE_Sensitivity_noapop.mat")
figure;
bar(categorical(regexprep(par_names(order),"_"," "),regexprep(par_names(order),"_"," ")),mu_star)