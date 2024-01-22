clearvars;

sm.fn = @computeTimeSeries;
model_type = ["exponential","logistic","von_bertalanffy"];
% model_type = "exponential";

save_fig_opts.save_figs = false;
save_fig_opts.file_types = ["fig","png"];



addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../SurrogateModelFns/")

nsamps = 3;
files.data = "../PostAnalysis/data/summary.mat";
for i = 1:length(model_type)
    load(sprintf("data/OptimalParameters_%s.mat",model_type(i)),"fstar")
    fstar = fstar(:);
    [~,order] = sort(fstar,"ascend");
    opts.abm_vec_inds(i) = order(1);
end
% 
% 
% % linear, exp growth, logistic
% opts.abm_vec_inds = [18,40,78];

for i = 1:length(model_type)
    switch model_type(i)
        case "exponential"
            opts.par_names = "\lambda";
        case "logistic"
            opts.par_names = ["r","K"];
        case "von_bertalanffy"
            opts.par_names = ["\alpha","\nu","\beta"];
    end
    files.optimal_parameters = sprintf("data/OptimalParameters_%s.mat",model_type(i));
    sm.opts.model_type = model_type(i);
    f{i} = testSMFitToABM(files, nsamps, sm, opts);
    fig_names_spec = ["SampleFitsOfSMToABM_%s","BestSMParameterDistributions_%s"];
    for j = numel(fig_names_spec):-1:1
        save_fig_opts.fig_names(j) = sprintf(fig_names_spec(j),model_type(i));
    end
end