function plotFunction_EndStep(M)

%% cell locations
if M.plot_pars.plot_location

    M.fig.scatter_plots(1).XData = M.healthy(:,M.I.subs(1))-M.grid.center(1);
    M.fig.scatter_plots(1).YData = M.healthy(:,M.I.subs(2))-M.grid.center(2);

    M.fig.scatter_plots(2).XData = M.tumor(:,M.I.subs(1))-M.grid.center(1);
    M.fig.scatter_plots(2).YData = M.tumor(:,M.I.subs(2))-M.grid.center(2);
    
%     view(M.fig.ax(M.fig.scatter_ind),[36*M.t 30 40])
end

%% tracked quantities
%% population plot
M.fig.population_plots(1).XData(end+1) = M.t;
M.fig.population_plots(1).YData(end+1) = M.healthy_count;
M.fig.population_plots(2).XData(end+1) = M.t;
M.fig.population_plots(2).YData(end+1) = M.tumor_count;

%% subpopulations proportion plot
M.fig.subpop_plots(1).XData(end+1) = M.t;
M.fig.subpop_plots(1).YData(end+1) = M.tracked.healthy_types(M.i,1)/M.tracked.healthy_count(M.i-1);
M.fig.subpop_plots(2).XData(end+1) = M.t;
M.fig.subpop_plots(2).YData(end+1) = M.tracked.healthy_types(M.i,2)/M.tracked.healthy_count(M.i-1);
M.fig.subpop_plots(3).XData(end+1) = M.t;
M.fig.subpop_plots(3).YData(end+1) = M.tracked.healthy_types(M.i,3)/M.tracked.healthy_count(M.i-1);
M.fig.subpop_plots(4).XData(end+1) = M.t;
M.fig.subpop_plots(4).YData(end+1) = M.tracked.tumor_types(M.i,1)/M.tracked.tumor_count(M.i-1);
M.fig.subpop_plots(5).XData(end+1) = M.t;
M.fig.subpop_plots(5).YData(end+1) = M.tracked.tumor_types(M.i,2)/M.tracked.tumor_count(M.i-1);
M.fig.subpop_plots(6).XData(end+1) = M.t;
M.fig.subpop_plots(6).YData(end+1) = M.tracked.tumor_types(M.i,3)/M.tracked.tumor_count(M.i-1);

%% event plots
M.fig.events_plots(1).XData(end+1) = M.t;
M.fig.events_plots(1).YData(end+1) = M.tracked.healthy_prolif(M.i) / (M.healthy_count*M.dt);
M.fig.events_plots(2).XData(end+1) = M.t;
M.fig.events_plots(2).YData(end+1) = M.tracked.healthy_contact_inhibition(M.i) / (M.healthy_count*M.dt);
M.fig.events_plots(3).XData(end+1) = M.t;
M.fig.events_plots(3).YData(end+1) = M.tracked.healthy_apop(M.i) / (M.healthy_count*M.dt);

M.fig.events_plots(4).XData(end+1) = M.t;
M.fig.events_plots(4).YData(end+1) = M.tracked.tumor_prolif(M.i) / (M.tumor_count*M.dt);
M.fig.events_plots(5).XData(end+1) = M.t;
M.fig.events_plots(5).YData(end+1) = M.tracked.tumor_contact_inhibition(M.i) / (M.tumor_count*M.dt);
M.fig.events_plots(6).XData(end+1) = M.t;
M.fig.events_plots(6).YData(end+1) = M.tracked.tumor_apop(M.i) / (M.tumor_count*M.dt);

drawnow