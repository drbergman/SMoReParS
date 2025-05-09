% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

sm.fn = @computeTimeSeries;
% model_type = "exponential";
% model_type = "logistic";
model_type = "von_bertalanffy";

sm.opts.model_type = model_type;

save_fig_opts.save_figs = false;
save_fig_opts.file_types = ["fig","png"];

fig_names_spec = ["SampleFitsOfSMToABM_%s","BestSMParameterDistributions_%s"];
for i = numel(fig_names_spec):-1:1
    save_fig_opts.fig_names(i) = sprintf(fig_names_spec(i),model_type);
end

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../../SurrogateModelFns/")

switch model_type
    case "exponential"
        par_names = "\lambda";
    case "logistic"
        par_names = ["r","K"];
    case "von_bertalanffy"
        par_names = ["\alpha","\nu","\beta"];
end

nsamps = 1;
files.optimal_parameters = sprintf("data/OptimalParameters_%s.mat",model_type);
files.data = "../PostAnalysis/data/summary.mat";

opts.par_names = par_names;
opts.abm_vec_inds = 1;
f = testSMFitToABM(files,nsamps,sm,opts);

% saveFigures(f,save_fig_opts)
