function M = initializeHealthy(M)

all_subs = allCombos(M.grid.xx,M.grid.yy);

M.healthy = zeros(M.setup.N0,M.I.proliferation_timer);
d = sqrt(sum((all_subs-M.grid.center).^2,2));

% M.healthy(:,M.I.subs) = all_subs(round(linspace(1,M.V_tot,M.setup.N0)'),:);
M.healthy(:,M.I.subs) = all_subs(randperm(M.V_tot,M.setup.N0)',:);
% M.healthy(:,M.I.subs) = all_subs(datasample(1:M.V_tot,M.setup.N0,'Replace',false,'Weights',exp(M.setup.c*((d/radius).^M.setup.e))),:);

M.healthy(:,M.I.event) = 0; % event index

%% determining remaining wait times for proliferation on healthy cells (see my Onenote document in ABM>Rates--Exponential>Going backwards in time)
prolif_rate = M.healthy_pars.prolif_rate; % assume that all mutant cells have phiD value of gammaT at start for simplicity
u = rand(M.setup.N0,1);
i1 = u>=1/(1+M.healthy_pars.min_prolif_wait*prolif_rate);
i2 = ~i1;
x = zeros(M.setup.N0,1);
x(i1) = u(i1)*(1+M.healthy_pars.min_prolif_wait*prolif_rate)/prolif_rate - 1/prolif_rate - M.healthy_pars.min_prolif_wait;
x(i2) = log((1+M.healthy_pars.min_prolif_wait*prolif_rate)*u(i2))/prolif_rate - M.healthy_pars.min_prolif_wait;
M.healthy(:,M.I.proliferation_timer) = max(0,x+M.healthy_pars.min_prolif_wait); % time until next possible proliferation (days)

%% figure out where healthy belong on grid
% (healthy location relative to start of lattice)/(length of each step) =
% num of lattice "right" of start; add 1 because start has index 1
M.healthy(:,M.I.ind) = sub2ind(M.grid.size,M.healthy(:,M.I.subs(1)),M.healthy(:,M.I.subs(2))); % linear indices

M.healthy_count = size(M.healthy,1);