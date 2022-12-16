function pars = basePars()

pars.max_dt = 4 / 24; % number of days per step

pars.cell_width = 20; % in micrometers; cell is about 20micrometers in diameter
pars.min_prolif_wait = 9/24; % number of days all cells must wait at minimum between proliferations
pars.max_tumor_size = Inf; % if tumor exceeds this, stop simulation

%% neighbor parameters
pars.occmax_2d = 6; % below this threshold, a tumor/immune cell can divide; at and above, too many neighbors and so doesn't proliferate
pars.occmax_3d = 20; % below this threshold, a tumor/immune cell can divide; at and above, too many neighbors and so doesn't proliferate

%% tumor parameters
pars.prolif_rate = 3; % proliferation rate of tumor cells (per day); from Global/Step18/test_determine_apoptosis.m saw that proliferation rate (r4 there) should be about 1.5, but want the minimum time before proliferation to be 9 hours, so need to increase proliferation rate here to make the average still be 1.5 proliferations / day
pars.apop_rate = .05; % max death rate of tumor cells (per day)
pars.move_rate_microns = 1 / (1/24); % move rate in microns per day

pars.mitosis_duration = 4/24; % duration of G2/mitosis in hours

%% chemo pars
pars.chemo_death_rate = 1;
