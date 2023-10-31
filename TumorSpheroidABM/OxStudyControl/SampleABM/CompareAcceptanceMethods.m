clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.resolution = '-r1200';

method_1_file = "data/AcceptedParameters_New_all_profiles_resampled.mat";
method_2_file = "data/AcceptedParameters_New_specified_parameter_combinations_resampled_1_2_3.mat";

experimental_data_file = "../ODEFitting/data/ExperimentalData_New.mat";
name_suffix = "New";
%% load accepted parameters
M(1) = load(method_1_file,"accepted_parameters","cohort_name","files");
M(2) = load(method_2_file,"accepted_parameters","cohort_name","files");

assert(isequal(M(1).cohort_name,M(2).cohort_name))
assert(isequal(M(1).files,M(2).files))

C = load(sprintf("../../data/%s/output.mat",M(1).cohort_name));
A = load(M(1).files.abm_data,"D");
E = load(experimental_data_file,"D","t","C");

n_conditions = numel(E.C);
nsamps = 100; % number of samples to use if plotting individual simulations
rss_on_log_scale = true;

ylim_total_cell_trajectories = [0 1000];

%% colors
color = lines(4);
color = color([1,3],:);
%% select accepted and rejected abm parameter vectors
for i = 1:2
    accepted{i} = A.D(:,M(i).accepted_parameters);
    rejected{i} = A.D(:,~M(i).accepted_parameters);
end
%% unify the runs
for i = 1:2
    accepted{i} = arrayify(accepted{i},"A",1);
    rejected{i} = arrayify(rejected{i},"A",1);
end

%% get total cell counts
for i = 1:2
    accepted{i} = sum(accepted{i},2);
    rejected{i} = sum(rejected{i},2);
    n_accepted(i) = sum(M(i).accepted_parameters,"all");
    n_rejected(i) = sum(~M(i).accepted_parameters,"all");
end
%% accepted vectors
for i = 1:2
    % accepted_vectors{i} = listAcceptedVectors(M(i).accepted_parameters);
end
%% patch plot just accepted
f = figureOnRight("Name","CompareAcceptedTrajectorySummaries_" + name_suffix);
ax = gobjects(n_conditions,1);

patch_plot_opts.patchPlotCoordsOptions.omit_nan = true;
patch_plot_opts.patchPlotCoordsOptions.split_sd = true;
patch_plot_opts.min_val = 0;
patch_plot_opts.LineWidth = 0.5;

p = gobjects(n_conditions,1,2);
l = gobjects(n_conditions,1,2);
data_line = gobjects(n_conditions,1);
for ci = 1:n_conditions % condition index
    for tsi = 1:1 % time series index
        for mi = 2:-1:1
            ax(ci,tsi) = subplot(n_conditions,1,r2c(n_conditions,1,[ci,tsi]));
            hold on;
            patch_plot_opts.Color = color(mi,:);
            % patch_plot_opts.DisplayName = "Accepted";
            [p(ci,tsi,1),l(ci,tsi,1)] = patchPlot(ax(ci,tsi),E.t,squeeze(accepted{mi}),patch_plot_opts);
            % patch_plot_opts.Color = colors(2,:);
            % patch_plot_opts.DisplayName = "Rejected";
            % [p(ci,tsi,2),l(ci,tsi,2)] = patchPlot(ax(ci,tsi),E.t,squeeze(rejected{i}),patch_plot_opts);
        end
        data_line(ci,tsi) = plot(E.t,E.D(ci).A,"Color","black","LineStyle","--","LineWidth",0.5,"Marker","*","MarkerSize",3,"DisplayName","Data");
    end
end
% legend_entries = [l(ishandle(l));data_line];
% L = legend(ax(1,1),legend_entries,"Location","best");
% L.Position = [0.7241    0.6300    0.2045    0.1655];
ylabel(ax(1,1),"Cell Count")
% title(ax(1,2),"G2/M Proportion")
% xlabel(ax(n_conditions,:),"Time (d)")
% if n_conditions>1
%     xticks(ax(1:(n_conditions-1),:),[])
% end
xlabel("Time (d)")
set(ax,"FontSize",8)
f.Units = "inches";
f.Position(3) = 2;
f.Position(4) = 1;
ylim(ylim_total_cell_trajectories)

%%
margin = struct("left",0.21,"right",.02,"top",.05,"bottom",0.3);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)

%% overlaid plots with low alpha
f = figureOnRight("Name","CompareSamplesByMethod_" + name_suffix);
ax = gca; hold on;

for i = 1:2
    k = min(nsamps,n_accepted(i));
    plot(ax,E.t,squeeze(accepted{i}(:,randperm(n_accepted(i),k))),"Color",[color(i,:),sqrt(1/k)],"LineWidth",0.5);
end
l_dat = plot(ax,E.t,E.D.A,"Color","black","LineStyle","--","LineWidth",0.5,"Marker","*","MarkerSize",3);

% legend(ax,[l_acc(1),l_rej(1),l_dat],["Accepted","Rejected","Data"]);
xlabel(ax,"Time (d)")
ylabel(ax,"Cell Count")
set(ax,"FontSize",8)
f.Units = "inches";
f.Position(3) = 2;
f.Position(4) = 1;
ylim(ylim_total_cell_trajectories)

%%
margin = struct("left",0.21,"right",.02,"top",.05,"bottom",0.3);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

saveFigures(f,save_fig_opts)

%% RSS breakdown
f = figureOnRight("Name","CompareRSSByMethod_" + name_suffix);
ax = gobjects(n_conditions,1);
RSS_accepted = cell(1,2);
for mi = 1:2
    RSS_accepted{mi} = sum(((accepted{mi} - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1);
end
for ci = 1:n_conditions % condition index
    ax(ci,1) = subplot(n_conditions,1,ci);
    hold on;
    if rss_on_log_scale
        [~,edges] = histcounts(log10(cat(3,RSS_accepted{1}(1,ci,:),RSS_accepted{2}(1,ci,:))));
        histogram(ax(ci,1),log10(RSS_accepted{1}(1,ci,:)),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Method 1","FaceColor",color(1,:))
        histogram(ax(ci,1),log10(RSS_accepted{2}(1,ci,:)),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Method 2","FaceColor",color(2,:))
    else
        [~,edges] = histcounts(cat(3,RSS_accepted{1}(1,ci,:),RSS_accepted{2}(1,ci,:))); %#ok<UNRCH>
        histogram(ax(ci,1),RSS_accepted{1}(1,ci,:),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Method 1","FaceColor",color(1,:))
        histogram(ax(ci,1),RSS_accepted{2}(1,ci,:),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Method 2","FaceColor",color(2,:))
        ax(ci,1).XLim(1) = 0;
    end
    ylabel(ax(ci,1),E.C{ci} + "\muM","FontWeight","bold")
end
normalizeXLims(ax)
ylabel("PDF")
if rss_on_log_scale
    xt = xticks(ax(1));
    set(ax,"XTickLabels",10.^xt)
end
% legend_entries = [reshape(l(1,2,:),[],1);data_line(1,2)];
% L = legend(ax(1,1)); %#ok<NASGU>
% L.Position = [0.7241    0.6300    0.2045    0.1655];
% title(ax(1,2),"G2/M Proportion")
xlabel(ax(n_conditions,:),"RSS")
set(ax,"FontSize",8)

f.Units = "inches";
f.Position(3) = 2;
f.Position(4) = 1;

%%
margin = struct("left",0.2,"right",.02,"top",.02,"bottom",0.29);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)

