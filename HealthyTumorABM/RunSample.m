clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%%
M = allBaseParameters();

%%
M.save_pars.dt = Inf; % set to Inf to not save anything; otherwise dt is in days

M.setup.grid_size_microns_x = 2000;
M.setup.grid_size_microns_y = 2000;

M.pars.max_dt = 2 / 24; % number of days per step

M.healthy_pars.min_prolif_wait = 12 / 24;
M.healthy_pars.prolif_rate = 1;
M.healthy_pars.apop_rate = 0.01;
M.healthy_pars.occmax = 2;

M.healthy_pars.min_prolif_wait = 8 / 24;
M.tumor_pars.prolif_rate = 2;
M.tumor_pars.apop_rate = 0.005;
M.tumor_pars.occmax = 7;

M.setup.censor_date = 100;
M.setup.N0 = 5e3;

M.plot_pars.plot_fig = true;
M.plot_pars.plot_location = true;
M.plot_pars.make_movie = true;

M = simPatient(M);
