% run the model on a grid of parameter values with nsamps_per_condition at
% each lattice point. to vary a parameter, set its value to a column
% vector. simCohort will look for all parameters with size("parameter",1)>1
% and vary these.

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%% cohort structure
cohort_pars.nsamps_per_condition = 6;
cohort_pars.min_parfor_num = 4e8;
cohort_pars.linkingFunction = @linkByName; % a function to link any parameters across the cohort (e.g. arrest coefficients or dosing parameters)

% parameters to link; each cell containts a string array of parameter values to link (if they are being varied); they must be varied over the same number of parameters
% cohort_pars.linkings = {["g1a_to_g1","g2a_to_g1"],... % recovery is identical for all arrested compartments
%                         ["arrest_ec50_g1","arrest_ec50_g2"]}; % ec50 for arrest is the same

cohort_pars.linkings = {["g1_to_s","s_to_g2","g2_to_m","m_to_g1"]};
cohort_pars.check_cohort_grab = true;
cohort_pars.previous_cohort_output_pattern = "data/cohort_*/output.mat";
cohort_pars.sim_function = @simPatient;
cohort_pars.n_per_batch = 3^3*6;
cohort_pars.update_timer_every = 8;
cohort_pars.parpool_options.resources = "Processes";

%%
M = allBaseParameters();
%%

% current picked values for varied parameters for a randomly accepted
% control ABM parameter at index I=2049 of the 3^7 grid

M.setup.ndims = 2;
M.setup.censor_date = 3;
M.setup.N0 = 1e2;
M.setup.agent_initialization_location = "uniform";
% M.setup.carrying_capacity = linspace(500,1500,7)';
M.setup.carrying_capacity = 500;
% M.setup.carrying_capacity = 1167.7;
M.setup.use_rates_for_intitial_proportions = false;

M.save_pars.make_save = true;
M.save_pars.save_constants = false;
M.save_pars.dt = Inf;
M.save_pars.interpolate_tracked = true; % only save certain interpolated values
M.save_pars.t_min = (0:1:72) * 60;
M.save_pars.fields_to_keep = ["t","phases"];

M.pars.max_dt = 0.25 / 24; % number of days per step
M.pars.occmax_3d = 20;
% M.pars.occmax_2d = [4;5;6];
M.pars.occmax_2d = 5;
M.pars.apop_rate = 0;
% M.pars.move_rate_microns = 10 * linspace(0,2,7)';
M.pars.move_rate_microns = 8.791732909379968;
% M.pars.move_rate_microns = 8.791732909379968;

M.flags.arrest_is_death = false;

M.chemo_pars.concentration = 0;

M.chemo_pars.apoptosis_function = "hill";
M.chemo_pars.apop_c0 = 0;
% M.chemo_pars.apop_c1 = [0.1;0.8;1.5];
% M.chemo_pars.apop_ec50 = [0.5;3.5;6.5];

transition_factors = linspace(0.5,2,7)';
% transition_factors = 1; % the model was the least sensitive to these
% 
M.cycle_pars.g1_to_s = 2.10273160861396490872721187770366668701171875;
M.cycle_pars.s_to_g2 = 2.939904610492844572178228190750814974308013916015625;
M.cycle_pars.g2_to_m = 5.96565977742448705356537175248377025127410888671875;
M.cycle_pars.m_to_g1 = 24.038155802861748355780946440063416957855224609375;

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = false;

%%
simCohort(M,cohort_pars);

