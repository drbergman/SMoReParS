% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = ["SampleFitsOfSMToABM_LMS","BestSMParameterDistributions_LMS","RSSOfSMFitsToABM_LMS"];

files.par_file = "data/SMFitToABM_LMS.mat";

load(files.par_file,"cohort_name")

files.data_file = sprintf("../../data/%s/summary.mat",cohort_name);
files.sm_fit_file = "data/SMFitToData_LMS.mat";

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ODEFittingFns/")

load("../ODEFitting/data/SMFitToData_LMS.mat","fixed_pars","model_type");
D = parameterDisplayNameDictionary(model_type);
sm_par_display_names = ["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"low_dose_apop";"delta_dose_apop";"rho0"];

column_names = {["Control Count";"NOT Fit"],["Control G2/M Prop";"NOT Fit"];
                "0.75\mu M Count","0.75\mu M G2/M Prop";
                "7.55\mu M Count","7.55\mu M G2/M Prop"};

nsamps = 10;

load("data/SMFitToData_LMS.mat","fn","fn_opts")

%% set up parameter names
[~,I] = setdiff(sm_par_display_names,fixed_pars);
sm_par_display_names = sm_par_display_names(sort(I));
for i = 1:numel(sm_par_display_names)
    sm_par_display_names(i) = D(sm_par_display_names(i));
end

%% Test the fitting
[f,I] = testSMFitToABM(files,nsamps,fn,fn_opts,sm_par_display_names,column_names);

%% save the figures
saveFigures(f,save_fig_opts)

%% reset path
rmpath("../../../ODEFittingFns/")


