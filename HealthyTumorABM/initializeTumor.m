function M = initializeTumor(M)

if M.setup.set_tumor_time==0
    all_subs = allCombos(M.grid.xx,M.grid.yy);

    M.tumor = zeros(M.setup.N0,M.I.proliferation_timer);
    d = sqrt(sum((all_subs-M.grid.center).^2,2));

    % M.tumor(:,M.I.subs) = all_subs(round(linspace(1,M.V_tot,M.setup.N0)'),:);
    M.tumor(:,M.I.subs) = all_subs(randperm(M.V_tot,M.setup.N0)',:);
    % M.tumor(:,M.I.subs) = all_subs(datasample(1:M.V_tot,M.setup.N0,'Replace',false,'Weights',exp(M.setup.c*((d/radius).^M.setup.e))),:);

    M.tumor(:,M.I.event) = 0; % event index

    %% determining remaining wait times for proliferation on tumor cells (see my Onenote document in ABM>Rates--Exponential>Going backwards in time)
    prolif_rate = M.tumor_pars.prolif_rate; % assume that all mutant cells have phiD value of gammaT at start for simplicity
    u = rand(M.setup.N0,1);
    i1 = u>=1/(1+M.tumor_pars.min_prolif_wait*prolif_rate);
    i2 = ~i1;
    x = zeros(M.setup.N0,1);
    x(i1) = u(i1)*(1+M.tumor_pars.min_prolif_wait*prolif_rate)/prolif_rate - 1/prolif_rate - M.tumor_pars.min_prolif_wait;
    x(i2) = log((1+M.tumor_pars.min_prolif_wait*prolif_rate)*u(i2))/prolif_rate - M.tumor_pars.min_prolif_wait;
    M.tumor(:,M.I.proliferation_timer) = max(0,x+M.tumor_pars.min_prolif_wait); % time until next possible proliferation (days)

    %% figure out where tumor belong on grid
    % (tumor location relative to start of lattice)/(length of each step) =
    % num of lattice "right" of start; add 1 because start has index 1
    M.tumor(:,M.I.ind) = sub2ind(M.grid.size,M.tumor(:,M.I.subs(1)),M.tumor(:,M.I.subs(2))); % linear indices


elseif M.t == 0
    M.tumor = zeros(0,M.I.proliferation_timer);

else % find all cells within radius of center and convert them to tumor cells
    healthy_subs = M.healthy(:,M.I.subs);
    healthy_distance = hypot(healthy_subs(:,1) - M.grid.center(1),healthy_subs(:,2) - M.grid.center(2));
    I = healthy_distance * M.pars.cell_width < M.setup.initial_tumor_radius_microns;
    M.tumor = M.healthy(I,:);
    M.healthy(I,:) = [];

    M.healthy_count = size(M.healthy,1);

end

M.tumor_count = size(M.tumor,1);
