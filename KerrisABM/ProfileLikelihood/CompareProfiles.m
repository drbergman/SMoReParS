clearvars;

addpath("../../ProfileLikelihoodFns/")

model_type = ["exponential","logistic","von_bertalanffy"];

save_fig_opts.save_figs = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = "SampleProfilesOfSMFromABM_" + model_type;

opts.abm_vec_inds = [51,57,35];
f = cell(length(model_type),1);
for i = 1:length(model_type)
    switch model_type(i)
        case "exponential"
            sm_par_display_names = "\lambda";
        case "logistic"
            sm_par_display_names = ["r","K"];
        case "von_bertalanffy"
            sm_par_display_names = ["\alpha","\nu","\beta"];
    end
    profile_file = sprintf("data/ProfileLikelihoods_%s.mat",model_type(i));
    nsamps = 3;
    [f{i},~] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);
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



rmpath("../../ProfileLikelihoodFns/")

