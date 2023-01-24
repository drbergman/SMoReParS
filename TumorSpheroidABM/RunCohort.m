clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%% cohort structure
cohort_pars.nsamps_per_condition = 6;
cohort_pars.min_parfor_num = 4;

%%
M = allBaseParameters();
%%

M.setup.ndims = 2;
M.setup.grid_size_microns_x = 2000;
M.setup.grid_size_microns_y = 2000;
M.setup.grid_size_microns_z = 2000;
M.setup.censor_date = 3;
M.setup.N0 = 1e2;
M.setup.agent_initialization_location = "uniform";
M.setup.carrying_capacity = [500;1000;1500];

M.save_pars.make_save = true;
M.save_pars.dt = Inf;

M.pars.max_dt = 0.25 / 24; % number of days per step
M.pars.occmax_3d = 20;
M.pars.occmax_2d = [4;5;6];
M.pars.apop_rate = [0.05;0.1;0.15];
M.pars.move_rate_microns = 10 * [0;1;2];

M.cycle_pars.g1_to_s = 24/11 * [0.8;1;1.2];
M.cycle_pars.s_to_g2 = 24/8 * [0.8;1;1.2];
M.cycle_pars.g2_to_m = 24/4 * [0.8;1;1.2];
M.cycle_pars.m_to_g1 = 24/1 * [0.8;1;1.2];

M.cycle_pars.dna_check_g1 = true;
M.cycle_pars.dna_check_s = false;
M.cycle_pars.dna_check_g2 = true;
M.cycle_pars.dna_check_m = false;

M.cycle_pars.arrest_prob_g1 = 0.00;
M.cycle_pars.arrest_prob_s = 0.00;
M.cycle_pars.arrest_prob_g2 = 0.00;
M.cycle_pars.arrest_prob_m = 0.00;

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = false;

%%
simCohort(M,cohort_pars);

