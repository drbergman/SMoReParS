% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

fn = @computeTimeSeries;
fn_opts.model_type = "von_bertalanffy";

save_fig_opts.save_figs = false;
save_fig_opts.file_types = ["fig","png"];

fig_names_spec = ["SampleFitsOfSMToABM_%s","BestSMParameterDistributions_%s"];
for i = numel(fig_names_spec):-1:1
    save_fig_opts.fig_names(i) = sprintf(fig_names_spec(i),fn_opts.model_type);
end

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../ODEFittingFns/")


par_names = ["\alpha","\nu","\beta"];

nsamps = 36;
par_file = sprintf("data/OptimalParameters_%s.mat",fn_opts.model_type);
data_file = "../PostAnalysis/data/summary.mat";


f = testSMFitToABM(par_file,data_file,nsamps,fn,fn_opts,par_names);

saveFigures(f,save_fig_opts)
