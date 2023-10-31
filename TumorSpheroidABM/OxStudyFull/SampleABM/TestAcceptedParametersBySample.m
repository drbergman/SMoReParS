% A script to explore how effective SMoRe ParS was a selecting the best
% parameters

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.resolution = '-r1200';

name_suffix = "LMS_2";
% acceptance_method = "single_best";
acceptance_method = "all_profiles_resampled";

experimental_data_file = "../ODEFitting/data/ExperimentalData.mat";
load("data/AcceptedParameters_" + name_suffix + "_" + acceptance_method + ".mat","accepted_parameters","cohort_name","files")

E = load(experimental_data_file,"D","t","C");
nsamps = 10; % number of samples to use if plotting individual simulations
rss_on_log_scale = true;
%% load and process ids
load(sprintf("../../data/%s/output.mat",cohort_name),"ids","lattice_parameters","nsamps_per_condition")
for i = 1:numel(lattice_parameters)
    if iscell(lattice_parameters(i).path)
        continue
    end
    if all(lattice_parameters(i).path==["chemo_pars","concentration"])
        condition_dim = i;
        break;
    end
end
parameter_dims = setdiff(1:length(lattice_parameters),condition_dim);
sample_dim = length(lattice_parameters)+1;

ids = permute(ids,[condition_dim,parameter_dims,sample_dim]);

%% select accepted and rejected abm parameter vectors
n_pars_accepted = sum(accepted_parameters,"all");
n_conditions = length(E.D);
accepted = zeros([size(E.D(1).A),n_conditions,n_pars_accepted*nsamps_per_condition]);
ad_ind = 1*ones(n_conditions,1);
n_rejected = sum(~accepted_parameters,"all");
rejected = zeros([size(E.D(1).A),n_conditions,n_rejected*nsamps_per_condition]);
re_ind = 1*ones(n_conditions,1);
par_ind = cell(1,length(parameter_dims));
for i = 1:numel(ids)
    S = load(sprintf("../../data/sims/%s/output_final.mat",ids(i)),"tracked");
    [cond_ind,par_ind{:},~] = ind2sub(size(ids),i);
    total = sum(S.tracked.phases,2);
    state2_prop = sum(S.tracked.phases(:,[3,4,7,8]),2)./total;
    if accepted_parameters(par_ind{:})
        accepted(:,1,cond_ind,ad_ind(cond_ind)) = total;
        accepted(:,2,cond_ind,ad_ind(cond_ind)) = state2_prop;
        ad_ind(cond_ind) = ad_ind(cond_ind)+1;
    else
        rejected(:,1,cond_ind,re_ind(cond_ind)) = total;
        rejected(:,2,cond_ind,re_ind(cond_ind)) = state2_prop;
        re_ind(cond_ind) = re_ind(cond_ind)+1;
    end
end

%% select accepted and rejected abm parameter vectors
n_pars_accepted = sum(accepted_parameters,"all");
n_conditions = length(E.D);
time_series_all = zeros([size(E.D(1).A),numel(ids)]);
par_ind = cell(1,length(parameter_dims));
for i = 1:numel(ids)
    S = load(sprintf("../../data/sims/%s/output_final.mat",ids(i)),"tracked");
    [cond_ind,par_ind{:},sample_ind] = ind2sub(size(ids),i);
    total = sum(S.tracked.phases,2);
    state2_prop = sum(S.tracked.phases(:,[3,4,7,8]),2)./total;
    time_series_all(:,1,i) = total;
    time_series_all(:,2,i) = state2_prop;
end
time_series_all = reshape(time_series_all,6,2,n_conditions,[],nsamps_per_condition);

%% patch plot
f = figureOnRight("Name","AcceptedVsRejectedTrajectories_" + name_suffix);
ax = gobjects(3,2);

patch_plot_opts.patchPlotCoordsOptions.omit_nan = true;
patch_plot_opts.patchPlotCoordsOptions.split_sd = true;
patch_plot_opts.min_val = 0;
patch_plot_opts.LineWidth = 0.5;

colors = lines(2);
p = gobjects(3,2,2);
l = gobjects(3,2,2);
data_line = gobjects(3,2);
for ci = 1:3 % condition index
    for tsi = 1:2 % time series index
        ax(ci,tsi) = subplot(3,2,r2c(3,2,[ci,tsi]));
        hold on;
        if tsi==2
            patch_plot_opts.max_val = 1; % don't let proportion go beyond 1
        else
            patch_plot_opts.max_val = Inf; % count is unbounded
        end
        patch_plot_opts.Color = colors(1,:);
        patch_plot_opts.DisplayName = "Accepted";
        [p(ci,tsi,1),l(ci,tsi,1)] = patchPlot(ax(ci,tsi),E.t,squeeze(accepted(:,tsi,ci,:)),patch_plot_opts);
        patch_plot_opts.Color = colors(2,:);
        patch_plot_opts.DisplayName = "Rejected";
        [p(ci,tsi,2),l(ci,tsi,2)] = patchPlot(ax(ci,tsi),E.t,squeeze(rejected(:,tsi,ci,:)),patch_plot_opts);
        data_line(ci,tsi) = plot(E.t,E.D(ci).A(:,tsi),"Color","black","LineStyle","--","LineWidth",0.5,"Marker","*","DisplayName","Data","MarkerSize",3);
    end
    if ci>1
        ylabel(ax(ci,1),"C = " + E.C{ci} + "\muM","FontWeight","bold")
    else
        ylabel(ax(ci,1),"Control","FontWeight","bold")
    end
end
normalizeYLims(ax(:,1))
normalizeYLims(ax(:,2))
% legend_entries = [reshape(l(1,2,:),[],1);data_line(1,2)];
% L = legend(ax(1,2),legend_entries);
% L.Position = [0.7241    0.6300    0.2045    0.1655];
title(ax(1,1),"Cell Count")
title(ax(1,2),"G2/M Fraction")
xlabel(ax(3,:),"Time (d)")
xticks(ax(1:2,:),[])
set(ax,"FontSize",8)
f.Units = "inches";
f.Position(3) = 2.5;
f.Position(4) = 3.5;

%% margins for this figure
margin = struct("left",.15,"right",.02,"top",.04,"bottom",.09);
spacing = struct("horizontal",0.1,"vertical",0.07);
uniformAxisSpacing(ax,margin,spacing);

%% save this figure
saveFigures(f,save_fig_opts)

%% RSS breakdown
f = figureOnRight("Name","RSSBreakdownBySamples_" + name_suffix + "_" + acceptance_method);
ax = gobjects(3,2);
RSS_accepted = sum(((accepted - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1);
RSS_rejected = sum(((rejected - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1);
for ci = 1:3 % condition index
    for tsi = 1:2 % time series index
        ax(ci,tsi) = subplot(3,2,r2c(3,2,[ci,tsi]));
        hold on;
        if rss_on_log_scale
            [~,edges] = histcounts(log10(cat(4,RSS_accepted(1,tsi,ci,:),RSS_rejected(1,tsi,ci,:))));
            histogram(ax(ci,tsi),log10(RSS_accepted(1,tsi,ci,:)),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Accepted")
            histogram(ax(ci,tsi),log10(RSS_rejected(1,tsi,ci,:)),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Rejected")
        else
            [~,edges] = histcounts(cat(4,RSS_accepted(1,tsi,ci,:),RSS_rejected(1,tsi,ci,:))); %#ok<UNRCH>
            histogram(ax(ci,tsi),RSS_accepted(1,tsi,ci,:),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Accepted")
            histogram(ax(ci,tsi),RSS_rejected(1,tsi,ci,:),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Rejected")
            ax(ci,tsi).XLim(1) = 0;
        end
    end
    if ci>1
        ylabel(ax(ci,1),"C = " + E.C{ci} + "\muM","FontWeight","bold")
    else
        ylabel(ax(ci,1),"Control","FontWeight","bold")
    end
end
normalizeXLims(ax)
if rss_on_log_scale
    xt = xticks(ax(1));
    set(ax,"XTickLabels",10.^xt)
end
% legend_entries = [reshape(l(1,2,:),[],1);data_line(1,2)];
% L = legend(ax(1,2)); %#ok<NASGU>
% L.Position = [0.7241    0.6300    0.2045    0.1655];
title(ax(1,1),"Cell Count")
title(ax(1,2),"G2/M Fraction")
xlabel(ax(3,:),"RSS")
set(ax,"FontSize",8)
yticks(ax,[])
xticks(ax(1:end-1,:),[])

f.Units = "inches";
f.Position(3) = 2.5;
f.Position(4) = 2.5;

%% margins for this figure
margin = struct("left",.07,"right",.02,"top",.06,"bottom",.12);
spacing = struct("horizontal",0.05,"vertical",0.08);
uniformAxisSpacing(ax,margin,spacing);

%% delete unused axis
delete(ax(1,2).Children)
ax(1,2).Color = f.Color;
ax(1,2).XColor = "none";
ax(1,2).YColor = "none";
ax(1,2).YAxis.Label.Color = "black";

%% save this figure
saveFigures(f,save_fig_opts)

%% RSS breakdown by most/least likely using total RSS
f = figureOnRight("Name","RSSBreakdownBySamples_PlaneThreshold_" + name_suffix + "_" + acceptance_method);
ax = gobjects(3,2);
% RSS_accepted = sum(((accepted - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1);
% RSS_rejected = sum(((rejected - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1);
RSS_all = sum(((time_series_all - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1);
% RSS_all = cat(4,RSS_accepted,RSS_rejected);
RSS_all = squeeze(RSS_all);
RSS_all_mean = mean(RSS_all,4); % average across samples
RSS_total = sum(RSS_all_mean,1:2,"omitnan"); % get total RSS across all time series
% RSS_total = sum(RSS_all,[1,2,4],"omitnan")/6;
[~,order] = sort(RSS_total);
NN = 70;
for ci = 1:3 % condition index
    for tsi = 1:2 % time series index
        ax(ci,tsi) = subplot(3,2,r2c(3,2,[ci,tsi]));
        hold on;
        if rss_on_log_scale
            [~,edges] = histcounts(log10(reshape(RSS_all(tsi,ci,:,:),[],1)));
            histogram(ax(ci,tsi),reshape(log10(RSS_all(tsi,ci,order(1:NN),:)),[],1),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Most Likely")
            histogram(ax(ci,tsi),reshape(log10(RSS_all(tsi,ci,order(NN+1:end),:)),[],1),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Least Likely")
        else
            [~,edges] = histcounts(reshape(RSS_all(tsi,ci,:,:),[],1)); %#ok<UNRCH>
            histogram(ax(ci,tsi),reshape(RSS_all(tsi,ci,order(1:NN),:),[],1),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Most Likely")
            histogram(ax(ci,tsi),reshape(RSS_all(tsi,ci,order(NN+1:end),:),[],1),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Least Likely")
            ax(ci,tsi).XLim(1) = 0;
        end
    end
    ylabel(ax(ci,1),E.C{ci} + "\muM","FontWeight","bold")
end
normalizeXLims(ax)
if rss_on_log_scale
    xt = xticks(ax(1));
    set(ax,"XTickLabels",10.^xt)
end
% legend_entries = [reshape(l(1,2,:),[],1);data_line(1,2)];
L = legend(ax(1,2)); %#ok<NASGU>
% L.Position = [0.7241    0.6300    0.2045    0.1655];
title(ax(1,1),"Total Cell Count")
title(ax(1,2),"G2/M Proportion")
xlabel(ax(3,:),"RSS")
set(ax,"FontSize",20)

saveFigures(f,save_fig_opts)

%% patch plot
f = figureOnRight("Name","ComparisonOfAccepted_" + name_suffix);
ax = gobjects(3,2);

time_series_all_best = time_series_all(:,:,:,order(1:n_pars_accepted),:);
time_series_all_best = reshape(time_series_all_best,6,2,3,[]);
patch_plot_opts.patchPlotCoordsOptions.omit_nan = true;
patch_plot_opts.patchPlotCoordsOptions.split_sd = true;
patch_plot_opts.min_val = 0;

colors = lines(2);
p = gobjects(3,2,2);
l = gobjects(3,2,2);
data_line = gobjects(3,2);
for ci = 1:3 % condition index
    for tsi = 1:2 % time series index
        ax(ci,tsi) = subplot(3,2,r2c(3,2,[ci,tsi]));
        hold on;
        if tsi==2
            patch_plot_opts.max_val = 1; % don't let proportion go beyond 1
        else
            patch_plot_opts.max_val = Inf; % count is unbounded
        end
        patch_plot_opts.Color = colors(1,:);
        patch_plot_opts.DisplayName = "SMoRe ParS";
        [p(ci,tsi,1),l(ci,tsi,1)] = patchPlot(ax(ci,tsi),E.t,squeeze(accepted(:,tsi,ci,:)),patch_plot_opts);
        patch_plot_opts.Color = colors(2,:);
        patch_plot_opts.DisplayName = "RSS Ranking";
        [p(ci,tsi,2),l(ci,tsi,2)] = patchPlot(ax(ci,tsi),E.t,squeeze(time_series_all_best(:,tsi,ci,:)),patch_plot_opts);
        data_line(ci,tsi) = plot(E.t,E.D(ci).A(:,tsi),"Color","black","LineStyle","--","LineWidth",2,"Marker","*","DisplayName","Data");
    end
    ylabel(ax(ci,1),E.C{ci} + "\muM","FontWeight","bold")
end
legend_entries = [reshape(l(1,2,:),[],1);data_line(1,2)];
L = legend(ax(1,2),legend_entries);
L.Position = [0.7241    0.6300    0.2045    0.1655];
title(ax(1,1),"Total Cell Count")
title(ax(1,2),"G2/M Proportion")
xlabel(ax(3,:),"Time (d)")
xticks(ax(1:2,:),[])
set(ax,"FontSize",20)

saveFigures(f,save_fig_opts)

%% RSS total
norm_method = "pdf";
f = figureOnRight("Name","RSSTotalBySample_" + name_suffix + "_" + acceptance_method + "_" + norm_method);
ax = gca;
hold on;
RSS_accepted = sum(((accepted - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1:3,"omitnan");
RSS_rejected = sum(((rejected - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1:3,"omitnan");
[~,edges] = histcounts([RSS_accepted(:);RSS_rejected(:)]);
histogram(ax,RSS_accepted(:),edges,"Normalization",norm_method,"EdgeColor","none","DisplayName","Accepted")
histogram(ax,RSS_rejected(:),edges,"Normalization",norm_method,"EdgeColor","none","DisplayName","Rejected")
ax.XLim(1) = 0;
% L = legend(ax,"Location","best");
% L.Position = [0.7241    0.6300    0.2045    0.1655];
title(ax,"Total RSS")
xlabel(ax,"RSS")
set(ax,"FontSize",8)
ylabel(ax,upper(norm_method))

f.Units = "inches";
f.Position(3) = 2.5;
f.Position(4) = 1;

%% margins for this figure
margin = struct("left",.13,"right",.07,"top",.15,"bottom",.3);
spacing = struct("horizontal",0.1,"vertical",0.07);
uniformAxisSpacing(ax,margin,spacing);

%% save this figure
saveFigures(f,save_fig_opts)

%% overlaid plots with low alpha
ax = gobjects(3,2,2);
f = figureOnRight;
for ci = 1:3
    for tsi = 1:2
        ax(ci,tsi,1) = subplot(3,2,r2c(3,2,[ci,tsi]));
        hold on;
    end
end
f(2) = figureOnRight;
for ci = 1:3
    for tsi = 1:2
        ax(ci,tsi,2) = subplot(3,2,r2c(3,2,[ci,tsi]));
        hold on;
    end
end

n_samples_accepted = size(accepted,4);
n_rejected = size(rejected,4);

colors = lines(2);
for ci = 1:3 % condition index
    for tsi = 1:2 % time series index
        k = min(nsamps,n_samples_accepted);
        plot(ax(ci,tsi,1),E.t,squeeze(accepted(:,tsi,ci,randperm(n_samples_accepted,k))),"Color",[colors(1,:),sqrt(1/k)])
        k = min(nsamps,n_rejected);
        plot(ax(ci,tsi,2),E.t,squeeze(rejected(:,tsi,ci,randperm(n_rejected,k))),"Color",[colors(2,:),sqrt(1/k)])
        plot(ax(ci,tsi,1),E.t,E.D(ci).A(:,tsi),"Color","black","LineStyle","--","LineWidth",2,"Marker","*")
        plot(ax(ci,tsi,2),E.t,E.D(ci).A(:,tsi),"Color","black","LineStyle","--","LineWidth",2,"Marker","*")
    end
end

saveFigures(f,save_fig_opts)
