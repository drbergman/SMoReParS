function tumor_pars = baseTumorPars()

tumor_pars.min_prolif_wait = 9/24; % number of days healthy cells must wait at minimum between proliferations
tumor_pars.max_healthy_size = Inf; % if healthy exceeds this, stop simulation

%% neighbor parameters
tumor_pars.occmax = 7; % below this threshold, a healthy/immune cell can divide; at and above, too many neighbors and so doesn't proliferate

%% healthy parameters

tumor_pars.prolif_rate = 1.5/(1-1.5*9/24); % proliferation rate of healthy cells (per day); from Global/Step18/test_determine_apoptosis.m saw that proliferation rate (r4 there) should be about 1.5, but want the minimum time before proliferation to be 9 hours, so need to increase proliferation rate here to make the average still be 1.5 proliferations / day
tumor_pars.apop_rate = .05; % max death rate of healthy cells (per day)
