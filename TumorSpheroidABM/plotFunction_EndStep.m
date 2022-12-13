function plotFunction_EndStep(M)

%% cell locations
if M.plot_pars.plot_location
    m_ind = M.tumor(:,M.I.phase)==M.val.phase_m;
    g0_ind = M.tumor(:,M.I.phase)==M.val.phase_g0;
    g1_ind = M.tumor(:,M.I.phase)==M.val.phase_g1;

    if M.setup.ndims == 3
        M.fig.scatter_plots(M.val.phase_m).XData = M.tumor(m_ind,M.I.subs(1))-M.grid.center(1);
        M.fig.scatter_plots(M.val.phase_m).YData = M.tumor(m_ind,M.I.subs(2))-M.grid.center(2);
        M.fig.scatter_plots(M.val.phase_m).ZData = M.tumor(m_ind,M.I.subs(3))-M.grid.center(3);

        M.fig.scatter_plots(M.val.phase_g0).XData = M.tumor(g0_ind,M.I.subs(1))-M.grid.center(1);
        M.fig.scatter_plots(M.val.phase_g0).YData = M.tumor(g0_ind,M.I.subs(2))-M.grid.center(2);
        M.fig.scatter_plots(M.val.phase_g0).ZData = M.tumor(g0_ind,M.I.subs(3))-M.grid.center(3);

        M.fig.scatter_plots(M.val.phase_g1).XData = M.tumor(g1_ind,M.I.subs(1))-M.grid.center(1);
        M.fig.scatter_plots(M.val.phase_g1).YData = M.tumor(g1_ind,M.I.subs(2))-M.grid.center(2);
        M.fig.scatter_plots(M.val.phase_g1).ZData = M.tumor(g1_ind,M.I.subs(3))-M.grid.center(3);

        view(M.fig.ax(M.fig.scatter_ind),[36*M.t 30 40])
    else
        M.fig.scatter_plots(M.val.phase_m).XData = M.tumor(m_ind,M.I.subs(1))-M.grid.center(1);
        M.fig.scatter_plots(M.val.phase_m).YData = M.tumor(m_ind,M.I.subs(2))-M.grid.center(2);

        M.fig.scatter_plots(M.val.phase_g0).XData = M.tumor(g0_ind,M.I.subs(1))-M.grid.center(1);
        M.fig.scatter_plots(M.val.phase_g0).YData = M.tumor(g0_ind,M.I.subs(2))-M.grid.center(2);

        M.fig.scatter_plots(M.val.phase_g1).XData = M.tumor(g1_ind,M.I.subs(1))-M.grid.center(1);
        M.fig.scatter_plots(M.val.phase_g1).YData = M.tumor(g1_ind,M.I.subs(2))-M.grid.center(2);
    end
end

if M.setup.ndims == 3
    %% slice
    z_mid_log_tumor = M.tumor(:,M.I.subs(3))==M.grid.center(3);
    M.fig.cell_slice_plot(1).XData = M.tumor(z_mid_log_tumor,M.I.subs(1));
    M.fig.cell_slice_plot(1).YData = M.tumor(z_mid_log_tumor,M.I.subs(2));
    M.fig.cell_slice_plot(1).ZData = ones(length(M.fig.cell_slice_plot(1).XData),1);

    view(M.fig.ax(M.fig.cell_slice_ind),2);

    %% projections
    if M.NT>1
        [N,~,~] = histcounts2(M.tumor(:,M.I.subs(1)),M.tumor(:,M.I.subs(2)),(M.grid.xx(1)-.5):(M.grid.xx(end)+.5),(M.grid.yy(1)-.5):(M.grid.yy(end)+.5),'Normalization','countdensity');
        if all(size(N)>1)
            M.fig.tum_density_plot.ZData = N';
        else
            M.fig.tum_density_plot.ZData = [];
        end
    end

end

%% tracked quantities
%% population plot
M.fig.population_plots(1).XData(end+1) = M.t;
M.fig.population_plots(1).YData(end+1) = M.NT;

%% subpopulations proportion plot
M.fig.subpop_plots(1).XData(end+1) = M.t;
M.fig.subpop_plots(1).YData(end+1) = M.tracked.phase_cell_days(M.i,1)/(M.dt*M.tracked.NT(M.i-1));
M.fig.subpop_plots(2).XData(end+1) = M.t;
M.fig.subpop_plots(2).YData(end+1) = M.tracked.phase_cell_days(M.i,2)/(M.dt*M.tracked.NT(M.i-1));
M.fig.subpop_plots(3).XData(end+1) = M.t;
M.fig.subpop_plots(3).YData(end+1) = M.tracked.phase_cell_days(M.i,3)/(M.dt*M.tracked.NT(M.i-1));

M.fig.subpop_plots(4).XData(end+1) = M.t;
M.fig.subpop_plots(4).YData(end+1) = M.tracked.simple_types(M.i,1)/M.tracked.NT(M.i-1);
M.fig.subpop_plots(5).XData(end+1) = M.t;
M.fig.subpop_plots(5).YData(end+1) = M.tracked.simple_types(M.i,2)/M.tracked.NT(M.i-1);

%% event plots
M.fig.events_plots(1).XData(end+1) = M.t;
M.fig.events_plots(1).YData(end+1) = M.tracked.tum_prolif(M.i) / (M.tracked.NT(M.i-1)*M.dt);
M.fig.events_plots(2).XData(end+1) = M.t;
M.fig.events_plots(2).YData(end+1) = M.tracked.tum_contact_inhibition(M.i) / (M.tracked.NT(M.i-1)*M.dt);
M.fig.events_plots(3).XData(end+1) = M.t;
M.fig.events_plots(3).YData(end+1) = M.tracked.tum_apop(M.i) / (M.tracked.NT(M.i-1)*M.dt);
M.fig.events_plots(4).XData(end+1) = M.t;
M.fig.events_plots(4).YData(end+1) = M.tracked.tum_chemo_death(M.i) / (M.tracked.NT(M.i-1)*M.dt);

drawnow