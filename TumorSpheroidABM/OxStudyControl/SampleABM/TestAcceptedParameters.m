% A script to explore how effective SMoRe ParS was a selecting the best
% parameters

clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.resolution = '-r1200';

name_suffix = "New";
% acceptance_method = "single_best";
acceptance_method = "all_profiles_resampled";
% acceptance_method = "specified_parameter_combinations_resampled";

experimental_data_file = "../ODEFitting/data/ExperimentalData_New.mat";
load("data/AcceptedParameters_" + name_suffix + "_" + acceptance_method + ".mat","accepted_parameters","cohort_name","files")
C = load(sprintf("../../data/%s/output.mat",cohort_name));
A = load(files.abm_data,"D");
E = load(experimental_data_file,"D","t","C");
n_conditions = numel(E.C);
nsamps = 100; % number of samples to use if plotting individual simulations
rss_on_log_scale = false;

ylim_total_cell_trajectories = [0 1000];

%% process abm parameter names
abm_par_names = setABMParameterNames(C.lattice_parameters);
n_abm_pars = length(abm_par_names);

%% colors
accept_color = [0    0.4470    0.7410];
reject_color = [0.8500    0.3250    0.0980];

%% select accepted and rejected abm parameter vectors
accepted = A.D(:,accepted_parameters);
rejected = A.D(:,~accepted_parameters);

%% unify the runs
accepted = arrayify(accepted,"A",1);
rejected = arrayify(rejected,"A",1);

%% get total cell counts
accepted = sum(accepted,2);
rejected = sum(rejected,2);
n_accepted = sum(accepted_parameters,"all");
n_rejected = sum(~accepted_parameters,"all");

%% accepted vectors
accepted_vectors = listAcceptedVectors(accepted_parameters);

%% all correlations
f = figure("Name","PairCorrelations_" + name_suffix + "_" + acceptance_method);
ax = gobjects(n_abm_pars,n_abm_pars);
[X,Y] = ndgrid(1:3,1:3);
for ri = 1:n_abm_pars
    for ci = 1:n_abm_pars
        if ri==ci
            continue;
        end
        ax(ri,ci) = subplot(n_abm_pars,n_abm_pars,r2c(n_abm_pars,n_abm_pars,[ri,ci]));
        xx = accepted_vectors(:,ci);
        yy = accepted_vectors(:,ri);
        temp = zeros(3);
        for ii = 1:3
            for jj = 1:3
                temp(ii,jj) = sum(xx==ii & yy==jj);
            end
        end
        scatter(ax(ri,ci),X,Y,temp,'k','filled');
    end
end


%% bar chart of accepted
f = figure("Name","ProportionAccepted_" + name_suffix + "_" + acceptance_method);
ax = gca;
cat_names = ["Accepted","Rejected"];
cc = categorical(cat_names);
bar(cc,[n_accepted,n_rejected]/numel(accepted_parameters));
ylabel("Proportion")
ax.FontSize = 20;

saveFigures(f,save_fig_opts)

%% bar chart of accepted by parameter at low-medium-high
f = figure("Name","ProportionAcceptedByValue_" + name_suffix + "_" + acceptance_method);
ax = gca;
cat_names = ["Low","Medium","High"];
cc = categorical(cat_names,cat_names);
admission_probabilities_by_value = computeAdmissionProbabilities(accepted_parameters);
bar(cc,admission_probabilities_by_value);
ylabel("Proportion")
ax.FontSize = 20;

saveFigures(f,save_fig_opts)

%% bar chart of accepted by low-medium-high
f = figure("Name","ProportionAcceptedByParameter_" + name_suffix + "_" + acceptance_method);
ax = gca;
% cat_names = {{'Carrying','Capacity'},"Contact Inhibition","Migration Rate","G1$\to$S Transition Rate","S$\to$G2 Transition Rate","G2$\to$M Transition Rate","M$\to$G1 Transition Rate"};
% cc = categorical(cat_names,cat_names);
val_colors = [135, 206, 235;65, 105, 225;70, 130, 180]/255;
admission_probabilities_by_value = computeAdmissionProbabilities(accepted_parameters);
b=bar(admission_probabilities_by_value);
for i = 1:3
    b(i).FaceColor = val_colors(i,:);
end
ylabel("Proportion")
ax.FontSize = 8;
xticks(ax,[]);
% ax.XAxis.TickLabelInterpreter = "latex";
% xticklabels(ax,cat_names)

f.Units = "inches";
f.Position(3) = 4.5;
f.Position(4) = 1;

%%
margin = struct("left",0.08,"right",.02,"top",.02,"bottom",0.05);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)

%% bubble chart
f = figure("Name","ProportionAcceptedByValueBubble_" + name_suffix + "_" + acceptance_method);
ax = gca;
ytick_names = ["Low","Medium","High"];
xtick_names = abm_par_names;
admission_probabilities_by_value = sqrt(computeAdmissionProbabilities(accepted_parameters));
for val_ind = 1:length(ytick_names)
    for par_ind = 1:length(xtick_names)
        viscircles([par_ind,val_ind],admission_probabilities_by_value(par_ind,val_ind))
    end
end
ax.FontSize = 20;

saveFigures(f,save_fig_opts)

%% patch plot just accepted
f = figureOnRight("Name","AcceptedTrajectorySummaries_" + name_suffix + "_" + acceptance_method);
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
        ax(ci,tsi) = subplot(n_conditions,1,r2c(n_conditions,1,[ci,tsi]));
        hold on;
        patch_plot_opts.Color = accept_color;
        patch_plot_opts.DisplayName = "Accepted";
        [p(ci,tsi,1),l(ci,tsi,1)] = patchPlot(ax(ci,tsi),E.t,squeeze(accepted),patch_plot_opts);
        % patch_plot_opts.Color = colors(2,:);
        % patch_plot_opts.DisplayName = "Rejected";
        % [p(ci,tsi,2),l(ci,tsi,2)] = patchPlot(ax(ci,tsi),E.t,squeeze(rejected),patch_plot_opts);
        data_line(ci,tsi) = plot(E.t,E.D(ci).A,"Color","black","LineStyle","--","LineWidth",0.5,"Marker","*","MarkerSize",3,"DisplayName","Data");
    end
end
legend_entries = [l(ishandle(l));data_line];
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
f.Position(3) = 1.5;
f.Position(4) = 1;
ylim(ylim_total_cell_trajectories)

margin = struct("left",0.3,"right",.04,"top",.05,"bottom",0.3);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

saveFigures(f,save_fig_opts)

%% patch plot
f = figureOnRight("Name","AcceptedVsRejectedTrajectorySummaries_" + name_suffix + "_" + acceptance_method);
ax = gobjects(n_conditions,1);

patch_plot_opts.patchPlotCoordsOptions.omit_nan = true;
patch_plot_opts.patchPlotCoordsOptions.split_sd = true;
patch_plot_opts.min_val = 0;

p = gobjects(n_conditions,1,2);
l = gobjects(n_conditions,1,2);
data_line = gobjects(n_conditions,1);
for ci = 1:n_conditions % condition index
    for tsi = 1:1 % time series index
        ax(ci,tsi) = subplot(n_conditions,1,r2c(n_conditions,1,[ci,tsi]));
        hold on;
        patch_plot_opts.Color = accept_color;
        patch_plot_opts.DisplayName = "Accepted";
        [p(ci,tsi,1),l(ci,tsi,1)] = patchPlot(ax(ci,tsi),E.t,squeeze(accepted),patch_plot_opts);
        patch_plot_opts.Color = reject_color;
        patch_plot_opts.DisplayName = "Rejected";
        [p(ci,tsi,2),l(ci,tsi,2)] = patchPlot(ax(ci,tsi),E.t,squeeze(rejected),patch_plot_opts);
        data_line(ci,tsi) = plot(E.t,E.D(ci).A,"Color","black","LineStyle","--","LineWidth",2,"Marker","*","DisplayName","Data");
    end
end
legend_entries = [reshape(l(ishandle(l)),[],1);data_line];
L = legend(ax(1,1),legend_entries,"Location","best");
% L.Position = [0.7241    0.6300    0.2045    0.1655];
ylabel(ax(1,1),"Total Cell Count")
% title(ax(1,2),"G2/M Proportion")
% xlabel(ax(n_conditions,:),"Time (d)")
% if n_conditions>1
%     xticks(ax(1:(n_conditions-1),:),[])
% end
xlabel("Time (d)")
set(ax,"FontSize",20)

saveFigures(f,save_fig_opts)

%% RSS breakdown
f = figureOnRight("Name","RSSBreakdown_" + name_suffix + "_" + acceptance_method);
ax = gobjects(n_conditions,1);
RSS_accepted = sum(((accepted - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1);
RSS_rejected = sum(((rejected - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1);
for ci = 1:n_conditions % condition index
    ax(ci,1) = subplot(n_conditions,1,ci);
    hold on;
    if rss_on_log_scale
        [~,edges] = histcounts(log10(cat(3,RSS_accepted(1,ci,:),RSS_rejected(1,ci,:))));
        histogram(ax(ci,1),log10(RSS_accepted(1,ci,:)),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Accepted","FaceColor",accept_color)
        histogram(ax(ci,tsi),log10(RSS_rejected(1,ci,:)),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Rejected","FaceColor",reject_color)
    else
        [~,edges] = histcounts(cat(3,RSS_accepted(1,ci,:),RSS_rejected(1,ci,:))); %#ok<UNRCH>
        histogram(ax(ci,1),RSS_accepted(1,ci,:),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Accepted")
        histogram(ax(ci,1),RSS_rejected(1,ci,:),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Rejected")
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
f.Position(3) = 1.5;
f.Position(4) = 1;

%%
margin = struct("left",0.26,"right",.02,"top",.02,"bottom",0.29);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)

%% RSS total
norm_method = "cdf";
f = figureOnRight("Name","RSSTotal_" + name_suffix + "_" + acceptance_method + "_" + norm_method);
ax = gca;
hold on;
RSS_accepted = sum(((accepted - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1:2,"omitnan");
RSS_rejected = sum(((rejected - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1:2,"omitnan");
[~,edges] = histcounts([RSS_accepted(:);RSS_rejected(:)]);
histogram(ax,RSS_accepted(:),edges,"Normalization",norm_method,"EdgeColor","none","DisplayName","Accepted","FaceColor",accept_color)
histogram(ax,RSS_rejected(:),edges,"Normalization",norm_method,"EdgeColor","none","DisplayName","Rejected","FaceColor",reject_color)
ax.XLim(1) = 0;
L = legend(ax,"Location","best");
% L.Position = [0.7241    0.6300    0.2045    0.1655];
title(ax,"Total RSS")
xlabel(ax,"RSS")
set(ax,"FontSize",20)
ylabel(ax,upper(norm_method))

saveFigures(f,save_fig_opts)

%% overlaid plots with low alpha
f = figureOnRight("Name","AcceptedSamples_" + name_suffix + "_" + acceptance_method);
ax = gca; hold on;
f(2) = figureOnRight("Name","RejectedSamples_" + name_suffix + "_" + acceptance_method);
ax(2) = gca; hold on;

k = min(nsamps,n_accepted);
plot(ax(1),E.t,squeeze(accepted(:,randperm(n_accepted,k))),"Color",[accept_color,2*sqrt(1/k)])
k = min(nsamps,n_rejected);
plot(ax(2),E.t,squeeze(rejected(:,randperm(n_rejected,k))),"Color",[reject_color,2*sqrt(1/k)])
plot(ax(1),E.t,E.D.A,"Color","black","LineStyle","--","LineWidth",2,"Marker","*")
plot(ax(2),E.t,E.D.A,"Color","black","LineStyle","--","LineWidth",2,"Marker","*")

set(ax,"FontSize",20)
xlabel(ax,"Time (d)")
ylabel(ax,"Total Cell Count")
title(ax(1),"Accepted")
title(ax(2),"Rejected")

normalizeYLims(ax)
saveFigures(f,save_fig_opts)

%% overlaid plots with low alpha
f = figureOnRight("Name","SamplesByAcceptedRejected_" + name_suffix + "_" + acceptance_method);
ax = gca; hold on;

k = min(nsamps,n_accepted);
l_acc = plot(ax,E.t,squeeze(accepted(:,randperm(n_accepted,k))),"Color",[accept_color,sqrt(1/k)],"LineWidth",0.5);
k = min(nsamps,n_rejected);
l_rej = plot(ax,E.t,squeeze(rejected(:,randperm(n_rejected,k))),"Color",[reject_color,sqrt(1/k)],"LineWidth",0.5);
l_dat = plot(ax,E.t,E.D.A,"Color","black","LineStyle","--","LineWidth",0.5,"Marker","*","MarkerSize",3);

% legend(ax,[l_acc(1),l_rej(1),l_dat],["Accepted","Rejected","Data"]);
set(ax,"FontSize",20)
xlabel(ax,"Time (d)")
ylabel(ax,"Cell Count")
set(ax,"FontSize",8)
f.Units = "inches";
f.Position(3) = 1.5;
f.Position(4) = 1;
ylim(ylim_total_cell_trajectories)

margin = struct("left",0.3,"right",.04,"top",.05,"bottom",0.3);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

saveFigures(f,save_fig_opts)

%% overlaid plots with low alpha
f = figureOnRight("Name","SamplesByAcceptedRejectedt10_" + name_suffix + "_" + acceptance_method);
ax = gca; hold on;

k = min(nsamps,n_accepted);
l_acc = plot(ax,E.t,squeeze(accepted(:,randperm(n_accepted,k))),"Color",[accept_color,sqrt(1/k)],"LineWidth",0.5);
k = min(nsamps,n_rejected);
l_rej = plot(ax,E.t,squeeze(rejected(:,randperm(n_rejected,k))),"Color",[reject_color,sqrt(1/k)],"LineWidth",0.5);
l_dat = plot(ax,E.t,E.D.A,"Color","black","LineStyle","--","LineWidth",0.5,"Marker","*","MarkerSize",3);

% legend(ax,[l_acc(1),l_rej(1),l_dat],["Accepted","Rejected","Data"]);
set(ax,"FontSize",20)
xlabel(ax,"Time (d)")
ylabel(ax,"Cell Count")
set(ax,"FontSize",8)
f.Units = "inches";
f.Position(3) = 1.5;
f.Position(4) = 1;
xlim(10/24 * [.8 1.2])

%%
margin = struct("left",0.24,"right",.03,"top",.05,"bottom",0.3);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)

%% overlaid plots with low alpha
f = figureOnRight("Name","SamplesByAcceptedRejectedt24_" + name_suffix + "_" + acceptance_method);
ax = gca; hold on;

k = min(nsamps,n_accepted);
l_acc = plot(ax,E.t,squeeze(accepted(:,randperm(n_accepted,k))),"Color",[accept_color,sqrt(1/k)],"LineWidth",0.5);
k = min(nsamps,n_rejected);
l_rej = plot(ax,E.t,squeeze(rejected(:,randperm(n_rejected,k))),"Color",[reject_color,sqrt(1/k)],"LineWidth",0.5);
l_dat = plot(ax,E.t,E.D.A,"Color","black","LineStyle","--","LineWidth",0.5,"Marker","*","MarkerSize",3);

% legend(ax,[l_acc(1),l_rej(1),l_dat],["Accepted","Rejected","Data"]);
set(ax,"FontSize",20)
xlabel(ax,"Time (d)")
ylabel(ax,"Cell Count")
set(ax,"FontSize",8)
f.Units = "inches";
f.Position(3) = 1.5;
f.Position(4) = 1;
xlim(24/24 * [.8 1.2])

%%
margin = struct("left",0.24,"right",.06,"top",.05,"bottom",0.3);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)

%% Internal functions

function admission_probabilities_by_value = computeAdmissionProbabilities(accepted_parameters)

sz = size(accepted_parameters);
n_pars = length(sz);
n_vals = unique(sz); % number of vals for each parameter
assert(numel(n_vals)==1) % make sure that each parameter had the same number of values sampled

N = numel(accepted_parameters)/n_vals; % how many values are sampled with fixing a single ABM (scalar) parameter
admission_probabilities_by_value = zeros(n_pars,n_vals);

for i = 1:n_pars
    temp = accepted_parameters;
    temp = permute(temp,[i,setdiff(1:n_pars,i)]);
    temp = reshape(temp,n_vals,[]);
    admission_probabilities_by_value(i,:) = sum(temp,2)/N;
end

end


function abm_par_names = setABMParameterNames(LP)

abm_par_names = strings(length(LP),1);
for i = 1:length(abm_par_names)
    switch LP(i).path{end}
        case 'carrying_capacity'
            abm_par_names(i) = "Carrying Capacity";
        case 'occmax_2d'
            abm_par_names(i) = "Contact Inhibition";
        case 'move_rate_microns'
            abm_par_names(i) = "Migration Rate";
        case 'g1_to_s'
            abm_par_names(i) = "r_{1S}";
        case 's_to_g2'
            abm_par_names(i) = "r_{S2}";
        case 'g2_to_m'
            abm_par_names(i) = "r_{2M}";
        case 'm_to_g1'
            abm_par_names(i) = "r_{M1}";
    end
end
end

function accepted_vectors = listAcceptedVectors(accepted_parameters)
I = find(accepted_parameters);
X = cell(ndims(accepted_parameters),1);
[X{:}] = ind2sub(size(accepted_parameters),I);
accepted_vectors = [X{:}];



end