% A script to explore how effective SMoRe ParS was a selecting the best
% parameters

clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig","png"];

% admission_method = "single_best";
admission_method = "all_profiles_resampled";

experimental_data_file = "../ODEFitting/data/ExperimentalData.mat";
load("data/AdmittedParameters_" + admission_method + ".mat","admitted_parameters","cohort_name","files")
A = load(files.abm_data_file,"D");
E = load(experimental_data_file,"D","t","C");
nsamps = 1000; % number of samples to use if plotting individual simulations

%% select admitted and rejected abm parameter vectors
admitted = A.D(:,admitted_parameters);
rejected = A.D(:,~admitted_parameters);

%% unify the runs
admitted = arrayify(admitted,"A",1);
rejected = arrayify(rejected,"A",1);

%% add t=0 time point:
admitted = cat(1,repmat([100,0.1],[1,1,size(admitted,3:ndims(admitted))]),admitted);
rejected = cat(1,repmat([100,0.1],[1,1,size(rejected,3:ndims(rejected))]),rejected);

%% patch plot
f = figureOnRight("Name","AdmittedVsRejectedTrajectorySummaries");
ax = gobjects(3,2);

patch_plot_opts.patchPlotCoordsOptions.omit_nan = true;
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
        patch_plot_opts.DisplayName = "Admitted";
        [p(ci,tsi,1),l(ci,tsi,1)] = patchPlot(ax(ci,tsi),E.t,squeeze(admitted(:,tsi,ci,:)),patch_plot_opts);
        patch_plot_opts.Color = colors(2,:);
        patch_plot_opts.DisplayName = "Rejected";
        [p(ci,tsi,2),l(ci,tsi,2)] = patchPlot(ax(ci,tsi),E.t,squeeze(rejected(:,tsi,ci,:)),patch_plot_opts);
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

%% RSS values
f = figureOnRight("Name","RSSBreakdown_" + admission_method);
ax = gobjects(3,2);
RSS_admitted = sum(((admitted - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1);
RSS_rejected = sum(((rejected - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1);
for ci = 1:3 % condition index
    for tsi = 1:2 % time series index
        ax(ci,tsi) = subplot(3,2,r2c(3,2,[ci,tsi]));
        hold on;
        [~,edges] = histcounts(cat(4,RSS_admitted(1,tsi,ci,:),RSS_rejected(1,tsi,ci,:)));
        histogram(ax(ci,tsi),RSS_admitted(1,tsi,ci,:),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Admitted")
        histogram(ax(ci,tsi),RSS_rejected(1,tsi,ci,:),edges,"Normalization","pdf","EdgeColor","none","DisplayName","Rejected")
        ax(ci,tsi).XLim(1) = 0;
    end
    ylabel(ax(ci,1),E.C{ci} + "\muM","FontWeight","bold")
end
% legend_entries = [reshape(l(1,2,:),[],1);data_line(1,2)];
L = legend(ax(1,2));
% L.Position = [0.7241    0.6300    0.2045    0.1655];
title(ax(1,1),"Total Cell Count")
title(ax(1,2),"G2/M Proportion")
xlabel(ax(3,:),"RSS")
set(ax,"FontSize",20)

saveFigures(f,save_fig_opts)

%% RSS total
norm_method = "pdf";
f = figureOnRight("Name","RSSTotal_" + admission_method + "_" + norm_method);
ax = gca;
hold on;
RSS_admitted = sum(((admitted - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1:3,"omitnan");
RSS_rejected = sum(((rejected - arrayify(E.D,"A",1))./arrayify(E.D,"S",1)).^2,1:3,"omitnan");
[~,edges] = histcounts([RSS_admitted(:);RSS_rejected(:)]);
histogram(ax,RSS_admitted(:),edges,"Normalization",norm_method,"EdgeColor","none","DisplayName","Admitted")
histogram(ax,RSS_rejected(:),edges,"Normalization",norm_method,"EdgeColor","none","DisplayName","Rejected")
ax.XLim(1) = 0;
L = legend(ax,"Location","best");
% L.Position = [0.7241    0.6300    0.2045    0.1655];
title(ax,"Total RSS")
xlabel(ax,"RSS")
set(ax,"FontSize",20)
ylabel(ax,upper(norm_method))

saveFigures(f,save_fig_opts)

%% overlaid plots with low alpha
ax = gobjects(3,2,2);
f(2) = figureOnRight;
for ci = 1:3
    for tsi = 1:2
        ax(ci,tsi,1) = subplot(3,2,r2c(3,2,[ci,tsi]));
        hold on;
    end
end
f(3) = figureOnRight;
for ci = 1:3
    for tsi = 1:2
        ax(ci,tsi,2) = subplot(3,2,r2c(3,2,[ci,tsi]));
        hold on;
    end
end

patch_plot_opts.patchPlotCoordsOptions.omit_nan = true;
patch_plot_opts.min_val = 0;

n_admitted = size(admitted,4);
n_rejected = size(rejected,4);

colors = lines(2);
for ci = 1:3 % condition index
    for tsi = 1:2 % time series index
        k = min(nsamps,n_admitted);
        plot(ax(ci,tsi,1),E.t,squeeze(admitted(:,tsi,ci,randperm(n_admitted,k))),"Color",[colors(1,:),sqrt(1/k)])
        k = min(nsamps,n_rejected);
        plot(ax(ci,tsi,2),E.t,squeeze(rejected(:,tsi,ci,randperm(n_rejected,k))),"Color",[colors(2,:),sqrt(1/k)])
        plot(ax(ci,tsi,1),E.t,E.D(ci).A(:,tsi),"Color","black","LineStyle","--","LineWidth",2,"Marker","*")
        plot(ax(ci,tsi,2),E.t,E.D(ci).A(:,tsi),"Color","black","LineStyle","--","LineWidth",2,"Marker","*")
    end
end
