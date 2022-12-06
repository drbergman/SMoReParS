clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%%
M = allBaseParameters();

%%
M.setup.ndims = 2;

M.save_pars.dt = Inf; % set to Inf to not save anything; otherwise dt is in days

M.setup.grid_size_microns_x = 1000;
M.setup.grid_size_microns_y = 1000;
M.setup.grid_size_microns_z = 1000;
M.pars.max_dt = 2 / 24; % number of days per step

M.pars.occmax = 23;
M.pars.prolif_rate = 2;
% M.pars.move_rate_microns = 0;

M.setup.censor_date = 20;
M.setup.N0 = 1e2;

M.plot_pars.plot_fig = true;
M.plot_pars.plot_location = true;

M = simPatient(M);
