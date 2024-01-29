clearvars;

sm.fn = @computeTimeSeries;
model_type = ["exponential","logistic","von_bertalanffy"];
% model_type = "exponential";

% for the first set of figures below
save_fig_opts.save_figs = false;
save_fig_opts.file_types = ["fig","png"];

% for the second set of figures below
save_compare_fig_opts.save_figs = true;
save_compare_fig_opts.file_types = ["fig","png"];
save_compare_fig_opts.resolution = "-r1200";



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

f = cell(length(model_type),1);
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

%%

g = gobjects(3,1);
ax = gobjects(3,1);
model_colors = lines(3);
[~,sort_order] = sort(opts.abm_vec_inds,"ascend");
model_types_sorted = flip(model_type(sort_order));
[~,inv_sort_order] = sort(sort_order);
inv_sort_order(sort_order) = 1:3;
for i = 1:3
    model_of_best = model_types_sorted(i);
    g(i) = figureOnRight("Name",sprintf("SMFitsOf_%s_best",model_of_best));
    ax(i) = gca;
    copyobj(f{1}(1).Children(i).Children(2:3),ax(i))
    for j = 1:3
        h_temp = copyobj(f{j}(1).Children(i).Children(1),ax(i));
        h_temp.Color = model_colors(j,:);
    end
    
end
xlim(ax,[0 75])

saveFigures(g,save_compare_fig_opts)

