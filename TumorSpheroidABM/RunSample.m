clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%%
M = allBaseParameters();

%%
M.setup.ndims = 2;
M.setup.censor_date = 3;
M.setup.N0 = 1e2;
M.setup.agent_initialization_location = "uniform";
M.setup.carrying_capacity = 1000;
M.setup.grid_size_microns_x = 2000;
M.setup.grid_size_microns_y = 2000;
M.setup.grid_size_microns_z = 2000;

M.save_pars.make_save = true;
M.save_pars.dt = Inf; % set to Inf to not save anything spatial; otherwise dt is in days
M.save_pars.interpolate_tracked = true; % only save certain interpolated values
M.save_pars.t_min = [0 10 24 36 48 72] * 60;
M.save_pars.fields_to_keep = ["t","phases"];


M.pars.max_dt = 0.25 / 24; % number of days per step
M.pars.occmax_3d = 20;
M.pars.occmax_2d = 5;
M.pars.move_rate_microns = 10;
M.pars.apop_rate = 0;

M.chemo_pars.concentration = 0.75;

M.chemo_pars.dna_check_g1 = true;
M.chemo_pars.dna_check_s = false;
M.chemo_pars.dna_check_g2 = true;
M.chemo_pars.dna_check_m = false;

M.chemo_pars.arrest_coeff_g1 = 0.01;
M.chemo_pars.arrest_coeff_s = 0.00;
M.chemo_pars.arrest_coeff_g2 = 0.00;
M.chemo_pars.arrest_coeff_m = 0.00;

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = true;
M.plot_pars.make_movie = false;

tic;
M = simPatient(M);
toc;

fprintf("Finished Simulation %s.\n",M.save_pars.sim_identifier)
