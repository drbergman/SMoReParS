clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%%
M = allBaseParameters();

%%
M.setup.ndims = 2;
M.setup.censor_date = 3;
M.setup.N0 = 1e2;
M.setup.agent_initialization_location = "uniform";
M.setup.carrying_capacity = 650;
M.setup.grid_size_microns_x = 2000;
M.setup.grid_size_microns_y = 2000;
M.setup.grid_size_microns_z = 2000;

M.save_pars.make_save = true;
M.save_pars.dt = 0.01; % set to Inf to not save anything; otherwise dt is in days

M.pars.max_dt = 0.25 / 24; % number of days per step
M.pars.chemo_death_rate = 0.05;
M.pars.occmax_3d = 20;
M.pars.occmax_2d = 5;
M.pars.move_rate_microns = 20;
M.pars.apop_rate = .05;

% M.cycle_pars.g1_to_s = 0.2;

M.cycle_pars.dna_check_g1 = false;
M.cycle_pars.dna_check_s = false;
M.cycle_pars.dna_check_g2 = false;
M.cycle_pars.dna_check_m = false;

M.cycle_pars.arrest_prob_g1 = 0.00;
M.cycle_pars.arrest_prob_s = 0.00;
M.cycle_pars.arrest_prob_g2 = 0.00;
M.cycle_pars.arrest_prob_m = 0.00;

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = true;
M.plot_pars.make_movie = false;

tic;
M = simPatient(M);
toc;
