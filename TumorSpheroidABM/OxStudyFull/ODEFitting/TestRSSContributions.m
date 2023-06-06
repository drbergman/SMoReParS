% This script will look at a particular fit of the SM to ABM output and
% break down the RSS of the n_conditions x n_time_series = 3
% x 2 = 6 contributions to the total RSS

clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];

load("data/SMFitToABM_FitAll.mat","P")
load("data/SMFitToData_FitAll.mat","fn_opts")
load("../../data/cohort_2305311216/summary.mat","D","C","n_conditions","n_time_series","t")

P = reshape(P,size(P,1),[]);

condition_titles = ["Control","0.75\mu M","7.55\mu M"];
time_series_titles = ["Count","G2/M Prop"];

RSS = zeros(n_conditions,n_time_series,size(P,2));
f = figureOnRight("Name","RSSBreakdown_SMFitToABM_FitAll");
ax = gobjects(n_conditions,n_time_series);
for j = 1:n_conditions
    title_str = condition_titles(j);
    for i = 1:size(P,2)
        sim_data = computeTimeSeries(P(:,i),t,C{j},fn_opts);
        RSS(j,:,i) = sum(((sim_data - D(j,i).A)./D(j,i).S).^2,1,"omitnan");
    end
    for k = 1:n_time_series
        ax(j,k) = subplot(n_conditions,n_time_series,r2c(n_conditions,n_time_series,[j,k]));
        histogram(ax(j,k),log10(RSS(j,k,:)))
        title(ax(j,k),condition_titles(j) + " " + time_series_titles(k))
    end
end
xlabel(ax(end,:),"RSS")
ylabel(ax(:,1),"Frequency")
normalizeXLims(f)
xt = xticks(ax(1));
set(ax,"XTickLabels",10.^xt)

saveFigures(f,save_fig_opts)
