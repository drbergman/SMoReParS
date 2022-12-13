function M = initializeFigure(M)

M.fig.handle = figure;
set(M.fig.handle,'units','normalized','outerposition',[0 0 1 1])

M.fig.nrows = 5;
M.fig.ncols = 6;

%% indices for the subplots in ax
M.fig.scatter_ind = 1;
M.fig.population_ind = 2;
M.fig.event_ind = 3;
M.fig.subpop_ind = 4;
M.fig.healthy_prob_ind = 5;
M.fig.tumor_prob_ind = 6;

ind_names = fieldnames(M.fig);
ind_names = ind_names(endsWith(ind_names,'ind'));

n_axes = 0;
for i = 1:length(ind_names)
    n_axes = max(n_axes,M.fig.(ind_names{i}));
end
M.fig.ax = gobjects(n_axes,1);

%% set up axes
cellfun(@(f) evalin('caller',[f ' = M.fig.' f ';']), fieldnames(M.fig));

scatter_locs = [1,3*ncols+4];

population_locs = 0*ncols + (5:6);
event_locs = 1*ncols + (5:6);
subpop_locs = 2*ncols + (5:6);

healthy_prob_locs = 4*ncols + (1:3);
tumor_prob_locs = 4*ncols + (4:6);


%% scatter plot
healthy_colors = winter(3);
tumor_colors = [0.9,0.1,0.1];
if M.plot_pars.plot_location
    M.fig.ax(scatter_ind) = subplot(nrows,ncols,scatter_locs);
    M.fig.ax(scatter_ind).Box = 'off';
    M.fig.ax(scatter_ind).NextPlot = 'add';
    M.fig.ax(scatter_ind).Color = "black";
    M.fig.ax(scatter_ind).XTick = [];
    M.fig.ax(scatter_ind).YTick = [];
    M.fig.ax(scatter_ind).ZTick = [];
    M.fig.ax(scatter_ind).XColor = 'none';
    M.fig.ax(scatter_ind).YColor = 'none';
    M.fig.ax(scatter_ind).ZColor = 'none';

    M.fig.scatter_plots(1) = scatter([],[],20,healthy_colors(1,:),'o','filled',...
        'DisplayName','Healthy Cells');
    M.fig.scatter_plots(2) = scatter([],[],20,tumor_colors(1,:),'o','filled',...
        'DisplayName','Tumor Cells');
    legend(M.fig.ax(M.fig.scatter_ind),'Location','NorthWest','AutoUpdate','off','Color',"white")

    M.fig.ax(M.fig.scatter_ind).Position(1) = 0.05;
    axis(M.fig.ax(M.fig.scatter_ind),[1,M.grid.size(1),1,M.grid.size(2)] - repelem(M.grid.center,1,2))

    axis square
%     M.fig.ax(M.fig.scatter_ind).Legend.Position(1) = M.fig.ax(M.fig.scatter_ind).Position(1) - M.fig.ax(M.fig.scatter_ind).Legend.Position(3);
end

%% time series plots
[xtick_vals,order] = sort(M.events.times);
symbols = {'I','^','.','C'};
vals = [0,1,2,Inf];
xx = [0 max(eps(),xtick_vals(1))];

[symbol_ind,~] = find((vals==M.events.event_index(order))');
xtick_labels = symbols(symbol_ind);

for i = length(xtick_vals):-1:2
    if xtick_vals(i)==xtick_vals(i-1)
        xtick_labels{i-1} = ';';
        xtick_labels(i) = [];
    end
end
xtick_vals = unique(xtick_vals);
   
%% population plot
M.fig.ax(population_ind) = subplot(nrows,ncols,population_locs);
M.fig.ax(population_ind).Title.String = 'Population Numbers';
M.fig.ax(population_ind).NextPlot = 'add';
M.fig.ax(population_ind).YLabel.String = 'Count';
M.fig.ax(population_ind).XAxis.TickValues = xtick_vals;
M.fig.ax(population_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(population_ind).XLim = xx;
M.fig.ax(population_ind).Color = [0 0 0];

M.fig.population_plots(1) = plot(M.t,M.setup.N0,'Color',healthy_colors(1,:),'LineStyle','-','LineWidth',2,'DisplayName','Healthy');
M.fig.population_plots(2) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','-','LineWidth',2,'DisplayName','Tumor');

legend(M.fig.ax(population_ind),'Location','west','AutoUpdate','off','Color',M.fig.handle.Color)

%% subpopulations proportion plot
M.fig.ax(subpop_ind) = subplot(nrows,ncols,subpop_locs);
M.fig.ax(subpop_ind).Title.String = 'Sub-Population Proportions';
M.fig.ax(subpop_ind).XLabel.String = 'Days';
M.fig.ax(subpop_ind).YLabel.String = 'Proportion';
M.fig.ax(subpop_ind).NextPlot = 'add';
% M.fig.ax(subpop_ind).XAxis.TickValues = xtick_vals;
% M.fig.ax(subpop_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(subpop_ind).XLim = xx;
% M.fig.ax(subpop_ind).YLim = [0 1];
M.fig.ax(subpop_ind).Color = [0 0 0];

M.fig.subpop_plots(1) = plot(M.t,0,'Color',healthy_colors(1,:),'LineStyle','-','LineWidth',2,'DisplayName','Healthy G0 Proportion');
M.fig.subpop_plots(2) = plot(M.t,0,'Color',healthy_colors(1,:),'LineStyle','--','LineWidth',2,'DisplayName','Healthy G1 Proportion');
M.fig.subpop_plots(3) = plot(M.t,0,'Color',healthy_colors(1,:),'LineStyle',':','LineWidth',2,'DisplayName','Healthy M Proportion');
M.fig.subpop_plots(4) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','-','LineWidth',2,'DisplayName','Tumor G0 Proportion');
M.fig.subpop_plots(5) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','--','LineWidth',2,'DisplayName','Tumor G1 Proportion');
M.fig.subpop_plots(6) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle',':','LineWidth',2,'DisplayName','Tumor M Proportion');

L = legend(M.fig.ax(subpop_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color);
L.Position(1) = 0.53;

%% event plot
M.fig.ax(event_ind) = subplot(nrows,ncols,event_locs);
M.fig.ax(event_ind).Title.String = 'Event Rates in Update';
M.fig.ax(event_ind).YLabel.String = 'days^{-1}';
M.fig.ax(event_ind).NextPlot = 'add';
M.fig.ax(event_ind).XAxis.TickValues = xtick_vals;
M.fig.ax(event_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(event_ind).XLim = xx;
M.fig.ax(event_ind).YScale = 'log';
M.fig.ax(event_ind).Color = [0 0 0];

M.fig.events_plots(1) = plot(M.t,0,'Color',healthy_colors(1,:),'LineStyle','-','LineWidth',2,'DisplayName','Healthy Proliferation');
M.fig.events_plots(2) = plot(M.t,0,'Color',healthy_colors(1,:),'LineStyle','--','LineWidth',2,'DisplayName','Healthy Contact Inhibitions');
M.fig.events_plots(3) = plot(M.t,0,'Color',healthy_colors(1,:),'LineStyle',':','LineWidth',2,'DisplayName','Healthy Apoptosis');
M.fig.events_plots(4) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','-','LineWidth',2,'DisplayName','Tumor Proliferation');
M.fig.events_plots(5) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','--','LineWidth',2,'DisplayName','Tumor Contact Inhibitions');
M.fig.events_plots(6) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle',':','LineWidth',2,'DisplayName','Tumor Apoptosis');

L = legend(M.fig.ax(event_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color);
L.Position(1) = 0.53;

%% not time series
%% healthy probabilities
M.fig.ax(healthy_prob_ind) = subplot(nrows,ncols,healthy_prob_locs);
max_healthy_prob = sum(1-exp(-M.healthy_pars.prolif_rate*M.pars.max_dt));
M.fig.ax(healthy_prob_ind).YLim = [0 max_healthy_prob];
M.fig.ax(healthy_prob_ind).Title.String = 'Healthy Outcome Probabilities';
M.fig.ax(healthy_prob_ind).NextPlot = 'replacechildren';

M.fig.healthy_prob_bar = bar(NaN,NaN(2,1),1,'stacked','EdgeAlpha',0);
M.fig.ax(M.fig.healthy_prob_ind).XTick = [];

events = {'prolif','apop'};
order = [2,1];
legend(M.fig.ax(healthy_prob_ind),flip(M.fig.healthy_prob_bar),flip(events(order)),'AutoUpdate','off','Location','northwest')

%% tumor probabilities
M.fig.ax(tumor_prob_ind) = subplot(nrows,ncols,tumor_prob_locs);
max_tumor_prob = sum(1-exp(-M.tumor_pars.prolif_rate*M.pars.max_dt));
M.fig.ax(tumor_prob_ind).YLim = [0 max_tumor_prob];
M.fig.ax(tumor_prob_ind).Title.String = 'Tumor Outcome Probabilities';
M.fig.ax(tumor_prob_ind).NextPlot = 'replacechildren';

M.fig.tumor_prob_bar = bar(NaN,NaN(2,1),1,'stacked','EdgeAlpha',0);
M.fig.ax(M.fig.tumor_prob_ind).XTick = [];

events = {'prolif','apop'};
order = [2,1];
legend(M.fig.ax(tumor_prob_ind),flip(M.fig.tumor_prob_bar),flip(events(order)),'AutoUpdate','off','Location','northwest')


drawnow