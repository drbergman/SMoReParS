% run the model on a grid of parameter values with nsamps_per_condition at
% each lattice point. to vary a parameter, set its value to a column
% vector. simCohort will look for all parameters with size("parameter",1)>1
% and vary these.

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%% cohort structure
cohort_pars.nsamps_per_condition = 6;
cohort_pars.min_parfor_num = 4;
cohort_pars.linkingFunction = @linkArrest; % a function to link any parameters across the cohort (e.g. arrest coefficients or dosing parameters)
cohort_pars.link_arrest_coeffs = false; % whether or not to even attempt to link arrest coefficients
cohort_pars.linkings = {}; % phases in a given cell have identical arrest coefficients (if they are being varied within the cohort)
cohort_pars.check_cohort_grab = true;
cohort_pars.previous_cohort_output_pattern = "data/cohort_*/output.mat";
cohort_pars.sim_function = @simPatient;
cohort_pars.update_timer_every = 10;
cohort_pars.parpool_options.resources = "Processes";

%%
M = allBaseParameters();
%%

M.setup.ndims = 2;
M.setup.censor_date = 3;
M.setup.N0 = 1e2;
M.setup.agent_initialization_location = "uniform";
M.setup.carrying_capacity = [500;1000;1500];
M.setup.use_rates_for_intitial_proportions = false;

M.save_pars.make_save = true;
M.save_pars.dt = Inf;
M.save_pars.interpolate_tracked = true; % only save certain interpolated values
M.save_pars.t_min = [0 10 24 36 48 72] * 60;
M.save_pars.fields_to_keep = ["t","phases"];

M.pars.max_dt = 0.25 / 24; % number of days per step
M.pars.occmax_3d = 20;
M.pars.occmax_2d = [4;5;6];
M.pars.apop_rate = 0.5;
M.pars.move_rate_microns = 10 * [0;1;2];

M.flags.apop_only_for_arrested = true;

% M.chemo_pars.concentration = [0;0.75;7.55];
M.chemo_pars.concentration = 7.55;
M.flags.arrest_is_death = false;

% transition_factors = [0.8;1;1.2];
transition_factors = 1;

M.cycle_pars.g1_to_s = 24/11 * transition_factors;
M.cycle_pars.s_to_g2 = 24/8 * transition_factors;
M.cycle_pars.g2_to_m = 24/4 * transition_factors;
M.cycle_pars.m_to_g1 = 24/1 * transition_factors;
M.cycle_pars.arrest_g1_to_g1 = 0.06;

M.chemo_pars.dna_check_g1 = true;
M.chemo_pars.dna_check_s = false;
M.chemo_pars.dna_check_g2 = true;
M.chemo_pars.dna_check_m = false;
M.chemo_pars.dna_check_g1a = false;

% arrest_coeffs = [0.025;0.05;0.075];

M.chemo_pars.arrest_function = "hill";

arrest_coeffs = 0.05;

M.chemo_pars.arrest_coeff_g1 = arrest_coeffs;
M.chemo_pars.arrest_coeff_s = 0;
M.chemo_pars.arrest_coeff_g2 = arrest_coeffs;
M.chemo_pars.arrest_coeff_m = 0;
M.chemo_pars.arrest_coeff_g1a = 0;

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = false;

%%
simCohort(M,cohort_pars);

