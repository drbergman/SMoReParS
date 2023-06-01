% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = ["SampleFitsOfSMToABM","BestSMParameterDistributions"];


cohort_name = "cohort_2305311216";

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ODEFittingFns/")

par_names = ["\alpha_R";"\alpha_P";"k_\alpha";"\delta_0";"k_\delta";"\rho_0"];


nsamps = 10;
par_file = "data/OptimalParameters.mat";
data_file = sprintf("../../data/%s/summary.mat",cohort_name);

fn = @computeTimeSeries;
model_type = "LogisticModel";
load("data/ODEFitToData.mat","fixed_pars")
D = parameterOrdering(model_type);

[p,~,~] = basePars(fixed_pars);
fixed_inds = zeros(numel(fixed_pars),1);
for i = 1:numel(fixed_pars)
    fixed_inds(i) = D(fixed_pars(i));
end
fixed_vals = p(fixed_inds);
p(fixed_inds) = [];
fn_opts.p_setup_fn = createParameterSetupFunction(model_type,fixed_pars);

[f,I] = testSMFitToABM(par_file,data_file,nsamps,fn,fn_opts,par_names);

saveFigures(f,save_fig_opts)

%% reset path
rmpath("../../../ODEFittingFns/")


