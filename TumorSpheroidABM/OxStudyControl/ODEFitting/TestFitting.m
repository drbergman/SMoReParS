% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = ["SampleFitsOfSMToABM","BestSMParameterDistributions"];

cohort_name = "cohort_230124175743017";

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ODEFittingFns/")

par_names = ["\lambda","\alpha","K"];

nsamps = 10;
par_file = "data/OptimalParameters_Using_DependentStatesAndTimeSeries.mat";
data_file = sprintf("../../data/%s/summary_short.mat",cohort_name);

fn = @computeTimeSeries;
fn_opts.condition_on_previous = false;

f = testSMFitToABM(par_file,data_file,nsamps,fn,fn_opts,par_names);

saveFigures(f,save_fig_opts)
