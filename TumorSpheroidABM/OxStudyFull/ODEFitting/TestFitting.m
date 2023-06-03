% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = ["SampleFitsOfSMToABM_Fit_b","BestSMParameterDistributions_Fit_b"];

cohort_name = "cohort_2305311216";

par_file = "data/SMFitToABM_Fit_b.mat";
data_file = sprintf("../../data/%s/summary.mat",cohort_name);

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ODEFittingFns/")

par_names = ["\alpha_R";"\alpha_P";"k_\alpha";"\delta_0";"k_\delta";"b";"\rho_0"];
column_names = {["Control Count";"NOT Fit"],["Control G2/M Prop";"NOT Fit"];
                "0.75\mu M Count","0.75\mu M G2/M Prop";
                "7.55\mu M Count","7.55\mu M G2/M Prop"};

nsamps = 10;

load("data/SMFitToData_Fit_b.mat","fixed_pars","fn","lb","ub","fn_opts","model_type")
[p,~,~,~] = fixParameters(model_type,fixed_pars);

[f,I] = testSMFitToABM(par_file,data_file,nsamps,fn,fn_opts,par_names,column_names);

saveFigures(f,save_fig_opts)

%% reset path
rmpath("../../../ODEFittingFns/")


