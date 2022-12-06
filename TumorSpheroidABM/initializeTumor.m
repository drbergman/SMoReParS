function M = initializeTumor(M)

if M.setup.ndims == 3
    radius = nthroot((3/4)*(1/pi())*M.setup.N0,3); % radius of sphere with volume equal to the number of cells
    all_subs = allCombos(M.grid.xx,M.grid.yy,M.grid.zz);

    M.tumor = zeros(M.setup.N0,M.I.proliferation_timer);
    d = sqrt(sum((all_subs-M.grid.center).^2,2));

    M.tumor(:,M.I.subs) = all_subs(datasample(1:M.V_tot,M.setup.N0,'Replace',false,'Weights',exp(M.setup.c*((d/radius).^M.setup.e))),:);

    M.tumor(:,M.I.event) = 0; % event index

    %% determining remaining wait times for proliferation on tumor cells (see my Onenote document in ABM>Rates--Exponential>Going backwards in time)
    prolif_rate = M.pars.prolif_rate; % assume that all mutant cells have phiD value of gammaT at start for simplicity
    u = rand(M.setup.N0,1);
    i1 = u>=1/(1+M.pars.min_prolif_wait*prolif_rate);
    i2 = ~i1;
    x = zeros(M.setup.N0,1);
    x(i1) = u(i1)*(1+M.pars.min_prolif_wait*prolif_rate)/prolif_rate - 1/prolif_rate - M.pars.min_prolif_wait;
    x(i2) = log((1+M.pars.min_prolif_wait*prolif_rate)*u(i2))/prolif_rate - M.pars.min_prolif_wait;
    M.tumor(:,M.I.proliferation_timer) = max(0,x+M.pars.min_prolif_wait); % time until next possible proliferation (days)

    %% figure out where tumors belong on grid
    % (tumor location relative to start of lattice)/(length of each step) =
    % num of lattice "right" of start; add 1 because start has index 1
    M.tumor(:,M.I.ind) = sub2ind(M.grid.size,M.tumor(:,M.I.subs(1)),M.tumor(:,M.I.subs(2)),M.tumor(:,M.I.subs(3))); % linear indices

else
    radius = nthroot((1/pi())*M.setup.N0,2); % radius of sphere with area equal to the number of cells
    all_subs = allCombos(M.grid.xx,M.grid.yy);

    M.tumor = zeros(M.setup.N0,M.I.proliferation_timer);
    d = sqrt(sum((all_subs-M.grid.center).^2,2));

    M.tumor(:,M.I.subs) = all_subs(datasample(1:M.V_tot,M.setup.N0,'Replace',false,'Weights',exp(M.setup.c*((d/radius).^M.setup.e))),:);

    M.tumor(:,M.I.event) = 0; % event index

    %% determining remaining wait times for proliferation on tumor cells (see my Onenote document in ABM>Rates--Exponential>Going backwards in time)
    prolif_rate = M.pars.prolif_rate; % assume that all mutant cells have phiD value of gammaT at start for simplicity
    u = rand(M.setup.N0,1);
    i1 = u>=1/(1+M.pars.min_prolif_wait*prolif_rate);
    i2 = ~i1;
    x = zeros(M.setup.N0,1);
    x(i1) = u(i1)*(1+M.pars.min_prolif_wait*prolif_rate)/prolif_rate - 1/prolif_rate - M.pars.min_prolif_wait;
    x(i2) = log((1+M.pars.min_prolif_wait*prolif_rate)*u(i2))/prolif_rate - M.pars.min_prolif_wait;
    M.tumor(:,M.I.proliferation_timer) = max(0,x+M.pars.min_prolif_wait); % time until next possible proliferation (days)

    %% figure out where tumors belong on grid
    % (tumor location relative to start of lattice)/(length of each step) =
    % num of lattice "right" of start; add 1 because start has index 1
    M.tumor(:,M.I.ind) = sub2ind(M.grid.size,M.tumor(:,M.I.subs(1)),M.tumor(:,M.I.subs(2))); % linear indices
end

M.NT = size(M.tumor,1);