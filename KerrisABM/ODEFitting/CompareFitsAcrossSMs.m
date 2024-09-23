clearvars;

sm.fn = @computeTimeSeries;
model_type = ["exponential","logistic","von_bertalanffy"];
% model_type = "exponential";

% for the first set of figures below
save_figs = true;
file_types = ["fig","png"];

resample_vb = true;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../SurrogateModelFns/")

nsamps = 3;
files.data = "../PostAnalysis/data/summary.mat";
RSS = cell(numel(model_type),1);
for i = 1:length(model_type)
    if model_type(i)~="von_bertalanffy" || ~resample_vb
        load(sprintf("data/OptimalParameters_%s.mat",model_type(i)),"fstar")
    else
        load(sprintf("data/OptimalParameters_%s_resampled.mat",model_type(i)),"fstar")
    end
    RSS{i} = fstar(:);
end

for i = 1:3
    if i~=3
        temp = RSS{i} - RSS{3};
    else
        temp = RSS{3} - RSS{2}; % compare vb to log
    end
    [~,opts.abm_vec_inds(i)] = min(temp);
    % [~,order] = sort(fstar,"ascend");
    % opts.abm_vec_inds(i) = order(1);
end
% opts.abm_vec_inds = [51; 66; 79]; % [best exp compared to vb; best log compared to vb; best vb compared to log] (when not down-sampling vb)
% opts.abm_vec_inds = [51; 57; 79]; % [best exp compared to vb; best log compared to vb; best vb compared to log] (when down-sampling vb)
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
    if model_type(i)~="von_bertalanffy" || ~resample_vb
        files.optimal_parameters = sprintf("data/OptimalParameters_%s.mat",model_type(i));
    else
        files.optimal_parameters = sprintf("data/OptimalParameters_%s_resampled.mat",model_type(i));
    end
    sm.opts.model_type = model_type(i);
    f{i} = testSMFitToABM(files, nsamps, sm, opts);
    fig_names_spec = ["SampleFitsOfSMToABM_%s","BestSMParameterDistributions_%s"];
    for j = numel(fig_names_spec):-1:1
        save_fig_opts.fig_names(j) = sprintf(fig_names_spec(j),model_type(i));
    end
end

%%

% for the second set of figures below
save_figs = true;
file_types = ["fig","png"];
resolution = "-r300";
reprint = true;

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
        h_temp.LineWidth = 0.5;
    end
    g(i).Units = "inches";
    g(i).Position(3:4) = [2,1.5];
end
xlim(ax,[0 75])
xticks(ax,0:25:75)
set(ax,"FontSize",8)
% set(ax,"YScale","log")
for i = 1:3
    if ax(i).YAxis.Scale == "log"
        g(i).Name = g(i).Name + "_log";
        ylim(ax(i),[1e2, 1e5])
    end
end

%%
saveFigures(g,save_figs=save_figs, file_types=file_types, resolution=resolution, reprint=reprint)

