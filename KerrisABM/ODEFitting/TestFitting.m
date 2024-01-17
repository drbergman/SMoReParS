% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

sm.fn = @computeTimeSeries;
% fn_opts.model_type = "von_bertalanffy";
sm.opts.model_type = "logistic";

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];

fig_names_spec = ["SampleFitsOfSMToABM_%s","BestSMParameterDistributions_%s"];
for i = numel(fig_names_spec):-1:1
    save_fig_opts.fig_names(i) = sprintf(fig_names_spec(i),fn_opts.model_type);
end

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../SurrogateModelFns/")

switch fn_opts.model_type
    case "logistic"
        par_names = ["r","K"];
    case "von_bertalanffy"
        par_names = ["\alpha","\nu","\beta"];
end

nsamps = 36;
par_file = sprintf("data/OptimalParameters_%s.mat",fn_opts.model_type);
data_file = "../PostAnalysis/data/summary.mat";


f = testSMFitToABM(par_file,data_file,nsamps,sm,par_names);

saveFigures(f,save_fig_opts)
