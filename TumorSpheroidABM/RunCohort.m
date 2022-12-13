clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%% cohort structure
cohort_pars.nsamps_per_condition = 6;
cohort_pars.min_parfor_num = 4;
cohort_pars.last_dose_is_no_dose = false;

%%
M = allBaseParameters();
%%

M.setup.ndims = 2;
M.setup.grid_size_microns_x = 2000;
M.setup.grid_size_microns_y = 2000;
M.setup.grid_size_microns_z = 2000;
M.setup.censor_date = 3;
M.setup.N0 = 1e3;
M.setup.agent_initialization_location = "uniform";
M.setup.carrying_capacity = [6e3;6500];

M.save_pars.dt = 0.25;

M.pars.max_dt = 4 / 24; % number of days per step
M.pars.chemo_death_rate = [0;0.2;0.5];
M.pars.occmax_3d = 20;
M.pars.occmax_2d = 6;
M.pars.min_prolif_wait = [9/24;12/24];
M.pars.prolif_rate = [1.8;2];
M.pars.apop_rate = [0.01;0.05];
M.pars.move_rate_microns = [0;20];

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = false;

%%
simCohort(M,cohort_pars);

% if cohort.nsamps_per_condition>=cohort.min_parfor_num
%     F(1:cohort.nsamps_per_condition) = parallel.FevalFuture;
%     ppool = gcp;
%     cohort.num_workers = ppool.NumWorkers;
%     for i = 1:cohort.nsamps_per_condition
%         F(i) = parfeval(ppool,@simPatient,1,M);
%     end
% else
%     cohort.num_workers = 1;
% end
% cohort.mu_n = 0;
% cohort.start = tic;
% cohort.batch_start = tic;
% 
% for si = cohort.nsamps_per_condition:-1:1
%     if cohort.nsamps_per_condition>=cohort.min_parfor_num
%         [idx,out_temp] = fetchNext(F);
%     else
%         idx = si;
%         out_temp = simPatient(M);
%     end
%     cohort = updateCohortTimer(cohort,cohort.nsamps_per_condition-si+1);
%     cohort.tracked(idx) = out_temp.tracked;
%     cohort.ids(idx) = out_temp.save_pars.sim_identifier;
% end
% 
% 
% save(sprintf("data/cohort_%d",tic),'-struct',"cohort")
