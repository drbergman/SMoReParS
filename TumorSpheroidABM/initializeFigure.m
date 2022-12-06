function M = initializeFigure(M)

M.fig.handle = figure;
set(M.fig.handle,'units','normalized','outerposition',[0 0 1 1])

M.fig.nrows = 5;
M.fig.ncols = 5;

%% indices for the subplots in ax
M.fig.scatter_ind = 1;
M.fig.tum_density_ind = 2;
M.fig.population_ind = 3;
M.fig.event_ind = 4;
M.fig.subpop_ind = 5;
M.fig.tum_prob_ind = 6;
M.fig.cell_slice_ind = 7;

ind_names = fieldnames(M.fig);
ind_names = ind_names(endsWith(ind_names,'ind'));

n_axes = 0;
for i = 1:length(ind_names)
    n_axes = max(n_axes,M.fig.(ind_names{i}));
end
M.fig.ax = gobjects(n_axes,1);

%% set up axes
cellfun(@(f) evalin('caller',[f ' = M.fig.' f ';']), fieldnames(M.fig));

scatter_locs = [1,2*ncols+3];
cell_slice_locs = [3*ncols+1,4*ncols+2];

tumor_density_locs = 0*ncols + 4;

population_locs = 1*ncols + (4:5);
event_locs = 3*ncols + (4:5);
subpop_locs = 2*ncols + (4:5);

tum_probs_locs = 0*ncols + 5;


%% scatter plot
tumor_colors = winter(3);
tumor_colormap = flipud(winter);

if M.plot_pars.plot_location
    M.fig.ax(scatter_ind) = subplot(nrows,ncols,scatter_locs);
    M.fig.ax(scatter_ind).Box = 'off';
    M.fig.ax(scatter_ind).NextPlot = 'add';
    M.fig.ax(scatter_ind).Color = M.fig.ax(scatter_ind).Parent.Color;
    M.fig.ax(scatter_ind).XTick = [];
    M.fig.ax(scatter_ind).YTick = [];
    M.fig.ax(scatter_ind).ZTick = [];
    M.fig.ax(scatter_ind).XColor = 'none';
    M.fig.ax(scatter_ind).YColor = 'none';
    M.fig.ax(scatter_ind).ZColor = 'none';

    M.fig.scatter_plots(1) = scatter3([],[],[],90,tumor_colors(1,:),'o','filled',...
        'DisplayName','Tumor Cells');
    view([0 30 40])
    legend(M.fig.ax(M.fig.scatter_ind),'Location','SouthWest','AutoUpdate','off','Color',"none")

    M.fig.ax(M.fig.scatter_ind).Legend.Position(1) = M.fig.ax(M.fig.scatter_ind).Position(1) - M.fig.ax(M.fig.scatter_ind).Legend.Position(3);
    axis(M.fig.ax(M.fig.scatter_ind),[1,M.grid.size(1),1,M.grid.size(2),1,M.grid.size(3)] - repelem(M.grid.center,1,2))

    axis square
end

%% slice plots
warning('off','MATLAB:contour:NonFiniteData')
%% cell slice plot
M.fig.ax(cell_slice_ind) = subplot(nrows,ncols,cell_slice_locs);
M.fig.ax(cell_slice_ind).NextPlot = 'add';
M.fig.ax(cell_slice_ind).XLim = [1,M.grid.size(1)];
M.fig.ax(cell_slice_ind).YLim = [1,M.grid.size(2)];
M.fig.cell_slice_plot(1) = scatter3([],[],[],60,tumor_colors(1,:),'o','filled',...
    'DisplayName','Tumor Cells','MarkerFaceAlpha',0.75);
M.fig.ax(cell_slice_ind).Title.String = 'Cell Slice (z=z_{mid})';
M.fig.ax(cell_slice_ind).XLabel.String = 'X coord (cells)';
M.fig.ax(cell_slice_ind).YLabel.String = 'Y coord (cells)';

axis square

if ~M.plot_pars.plot_location
    legend(M.fig.cell_slice_plot(1:6),'Location','bestoutside','AutoUpdate','off')
end

%% density projection plots
%% tumor density plot
M.fig.ax(tum_density_ind) = subplot(nrows,ncols,tumor_density_locs);
[~,M.fig.tum_density_plot] = contourf(M.grid.xx,M.grid.yy,NaN(M.grid.size(1:2))','LineColor','none');
M.fig.ax(tum_density_ind).Title.String = 'Tumor Projected onto (X,Y)';
M.fig.ax(tum_density_ind).XLabel.String = 'X coord (cells)';
M.fig.ax(tum_density_ind).YLabel.String = 'Y coord (cells)';
M.fig.tum_density_colorbar = colorbar;

colormap(M.fig.ax(tum_density_ind),tumor_colormap);

axis square

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
M.fig.ax(population_ind).XAxis.TickValues = xtick_vals;
M.fig.ax(population_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(population_ind).XLim = xx;
M.fig.ax(population_ind).Color = [0 0 0];

M.fig.population_plots(1) = plot(M.t,M.setup.N0,'w-','LineWidth',2,'DisplayName','N_T');
% M.fig.population_plots(2) = plot(M.t,M.setup.NI0,'b-','LineWidth',2,'DisplayName','N_I');

legend(M.fig.ax(population_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color)

%% subpopulations proportion plot
M.fig.ax(subpop_ind) = subplot(nrows,ncols,subpop_locs);
M.fig.ax(subpop_ind).Title.String = 'Sub-Population Proportions';
M.fig.ax(subpop_ind).NextPlot = 'add';
M.fig.ax(subpop_ind).XAxis.TickValues = xtick_vals;
M.fig.ax(subpop_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(subpop_ind).XLim = xx;
M.fig.ax(subpop_ind).YLim = [0 1];
M.fig.ax(subpop_ind).Color = [0 0 0];

M.fig.subpop_plots(1) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','-','LineWidth',2,'DisplayName','G0 Proportion');
M.fig.subpop_plots(2) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','--','LineWidth',2,'DisplayName','G1 Proportion');
M.fig.subpop_plots(3) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle',':','LineWidth',2,'DisplayName','M Proportion');

M.fig.subpop_plots(4) = plot(M.t,0,'Color',tumor_colors(2,:),'LineStyle','-','LineWidth',2,'DisplayName','Not proliferate');
M.fig.subpop_plots(5) = plot(M.t,0,'Color',tumor_colors(2,:),'LineStyle',':','LineWidth',2,'DisplayName','Proliferate');

legend(M.fig.ax(subpop_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color)

%% event plot
M.fig.ax(event_ind) = subplot(nrows,ncols,event_locs);
M.fig.ax(event_ind).Title.String = 'Event Rates in Update';
M.fig.ax(event_ind).XLabel.String = 'Days';
M.fig.ax(event_ind).NextPlot = 'add';
% M.fig.ax(event_ind).XAxis.TickValues = xtick_vals;
% M.fig.ax(event_ind).XAxis.TickLabels = xtick_labels;
M.fig.ax(event_ind).XLim = xx;
M.fig.ax(event_ind).YScale = 'log';
M.fig.ax(event_ind).Color = [0 0 0];

M.fig.events_plots(1) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','-','LineWidth',2,'DisplayName','Tumor Proliferation');
M.fig.events_plots(2) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle','--','LineWidth',2,'DisplayName','Tumor Contact Inhibitions');
M.fig.events_plots(3) = plot(M.t,0,'Color',tumor_colors(1,:),'LineStyle',':','LineWidth',2,'DisplayName','Tumor Apoptosis');

legend(M.fig.ax(event_ind),'Location','northwest','AutoUpdate','off','Color',M.fig.handle.Color)

%% not time series
%% tumor probabilities
M.fig.ax(tum_prob_ind) = subplot(nrows,ncols,tum_probs_locs);
max_tum_prob = sum(1-exp(-[M.pars.prolif_rate,M.pars.apop_rate,M.pars.move_rate]*M.pars.max_dt));
M.fig.ax(tum_prob_ind).YLim = [0 max_tum_prob];
M.fig.ax(tum_prob_ind).Title.String = 'Tum Outcome Probabilities';
M.fig.ax(tum_prob_ind).NextPlot = 'replacechildren';

M.fig.tum_prob_bar = bar(NaN,NaN(3,1),1,'stacked','EdgeAlpha',0);

events = {'prolif','apop','move'};
order = [3,2,1];
legend(M.fig.ax(tum_prob_ind),flip(M.fig.tum_prob_bar),flip(events(order)),'AutoUpdate','off','Location','northwest')

drawnow