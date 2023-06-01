% run the model on a grid of parameter values with nsamps_per_condition at
% each lattice point. to vary a parameter, set its value to a column
% vector. simCohort will look for all parameters with size("parameter",1)>1
% and vary these.

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%% cohort structure
cohort_pars.nsamps_per_condition = 6;
cohort_pars.min_parfor_num = 4;
cohort_pars.linkingFunction = @linkByName; % a function to link any parameters across the cohort (e.g. arrest coefficients or dosing parameters)

% parameters to link; each cell containts a string array of parameter values to link (if they are being varied); they must be varied over the same number of parameters
cohort_pars.linkings = {["g1a_to_g1","g2a_to_g1"],... % recovery is identical for all arrested compartments
                        ["arrest_ec50_g1","arrest_ec50_g2"]}; % ec50 for arrest is the same

cohort_pars.check_cohort_grab = true;
cohort_pars.previous_cohort_output_pattern = "data/cohort_*/output.mat";
cohort_pars.sim_function = @simPatient;
cohort_pars.n_per_batch = 3^4*6;
cohort_pars.update_timer_every = 8;
cohort_pars.parpool_options.resources = "Processes";

%%
M = allBaseParameters();
%%

M.setup.ndims = 2;
M.setup.censor_date = 3;
M.setup.N0 = 1e2;
M.setup.agent_initialization_location = "uniform";
% M.setup.carrying_capacity = [500;1000;1500];
M.setup.carrying_capacity = 1167.726550079491289579891599714756011962890625;
M.setup.use_rates_for_intitial_proportions = false;

M.save_pars.make_save = true;
M.save_pars.save_constants = false;
M.save_pars.dt = Inf;
M.save_pars.interpolate_tracked = true; % only save certain interpolated values
M.save_pars.t_min = [0 10 24 36 48 72] * 60;
M.save_pars.fields_to_keep = ["t","phases"];

M.pars.max_dt = 0.25 / 24; % number of days per step
M.pars.occmax_3d = 20;
% M.pars.occmax_2d = [4;5;6];
M.pars.occmax_2d = 5;
M.pars.apop_rate = 0;
% M.pars.move_rate_microns = 10 * [0;1;2];
M.pars.move_rate_microns = 8.791732909379968;

M.flags.arrest_is_death = false;

M.chemo_pars.concentration = [0;0.75;7.55];

M.chemo_pars.apoptosis_function = "hill";
M.chemo_pars.apop_c0 = 0;
M.chemo_pars.apop_c1 = [0.1;0.8;1.5];
M.chemo_pars.apop_ec50 = [0.5;3.5;6.5];

% % transition_factors = [0.8;1;1.2];
% transition_factors = 1; % the model was the least sensitive to these
% 
% M.cycle_pars.g1_to_s = 24/11 * transition_factors;
% M.cycle_pars.s_to_g2 = 24/8 * transition_factors;
% M.cycle_pars.g2_to_m = 24/4 * transition_factors;
% M.cycle_pars.m_to_g1 = 24/1 * transition_factors;

M.cycle_pars.g1_to_s = 2.10273160861396490872721187770366668701171875;
M.cycle_pars.s_to_g2 = 2.939904610492844572178228190750814974308013916015625;
M.cycle_pars.g2_to_m = 5.96565977742448705356537175248377025127410888671875;
M.cycle_pars.m_to_g1 = 24.038155802861748355780946440063416957855224609375;

recovery_rate = [0.03;0.06;0.09];
M.cycle_pars.g1a_to_g1 = recovery_rate;
M.cycle_pars.g2a_to_g1 = recovery_rate;

M.chemo_pars.dna_check_g1 = true;
M.chemo_pars.dna_check_g2 = true;

M.chemo_pars.arrest_function = "hill";
arrest_coeffs = [0.25;0.5;0.75];
arrest_ec50s = [0.5;3.5;6.5];

M.chemo_pars.arrest_coeff_g1 = arrest_coeffs;
M.chemo_pars.arrest_coeff_g2 = arrest_coeffs;
M.chemo_pars.arrest_ec50_g1 = arrest_ec50s;
M.chemo_pars.arrest_ec50_g2 = arrest_ec50s;

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = false;

%%
simCohort(M,cohort_pars);

