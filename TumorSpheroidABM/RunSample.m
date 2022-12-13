clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%%
M = allBaseParameters();

%%
M.setup.ndims = 2;

M.save_pars.dt = Inf; % set to Inf to not save anything; otherwise dt is in days

M.setup.grid_size_microns_x = 2000;
M.setup.grid_size_microns_y = 2000;
M.setup.grid_size_microns_z = 2000;
M.pars.max_dt = 4 / 24; % number of days per step

M.pars.chemo_death_rate = 1;

M.pars.occmax_3d = 20;
M.pars.occmax_2d = 6;
M.pars.prolif_rate = 2;
M.pars.move_rate_microns = 60;

M.setup.censor_date = 3;
M.setup.N0 = 1e3;
M.setup.agent_initialization_location = "uniform";

M.plot_pars.plot_fig = true;
M.plot_pars.plot_location = true;
M.plot_pars.make_movie = false;

tic;
M = simPatient(M);
toc;
