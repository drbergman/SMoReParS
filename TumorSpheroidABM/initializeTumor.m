function M = initializeTumor(M)

M.tumor = zeros(M.setup.N0,M.I.phase);
switch M.setup.agent_initialization_location
    case "center"
        if M.setup.ndims == 3
            radius = nthroot((3/4)*(1/pi())*M.setup.N0,3); % radius of sphere with volume equal to the number of cells
            all_subs = allCombos(M.grid.xx,M.grid.yy,M.grid.zz);
            d = sqrt(sum((all_subs-M.grid.center).^2,2));

            M.tumor(:,M.I.subs) = all_subs(datasample(1:M.V_tot,M.setup.N0,'Replace',false,'Weights',exp(M.setup.c*((d/radius).^M.setup.e))),:);

            %% figure out where tumors belong on grid
            % (tumor location relative to start of lattice)/(length of each step) =
            % num of lattice "right" of start; add 1 because start has index 1
            M.tumor(:,M.I.ind) = sub2ind(M.grid.size,M.tumor(:,M.I.subs(1)),M.tumor(:,M.I.subs(2)),M.tumor(:,M.I.subs(3))); % linear indices

        else
            radius = nthroot((1/pi())*M.setup.N0,2); % radius of sphere with area equal to the number of cells
            all_subs = allCombos(M.grid.xx,M.grid.yy);
            d = sqrt(sum((all_subs-M.grid.center).^2,2));

            M.tumor(:,M.I.subs) = all_subs(datasample(1:M.V_tot,M.setup.N0,'Replace',false,'Weights',exp(M.setup.c*((d/radius).^M.setup.e))),:);

            %% figure out where tumors belong on grid
            % (tumor location relative to start of lattice)/(length of each step) =
            % num of lattice "right" of start; add 1 because start has index 1
            M.tumor(:,M.I.ind) = sub2ind(M.grid.size,M.tumor(:,M.I.subs(1)),M.tumor(:,M.I.subs(2))); % linear indices
        end

    case "uniform"
        M.tumor(:,M.I.ind) = randperm(M.V_tot,M.setup.N0);
        if M.setup.ndims==3
            [M.tumor(:,M.I.subs(1)),M.tumor(:,M.I.subs(2)),M.tumor(:,M.I.subs(3))] = ind2sub(M.grid.size,M.tumor(:,M.I.ind));
        else
            [M.tumor(:,M.I.subs(1)),M.tumor(:,M.I.subs(2))] = ind2sub(M.grid.size,M.tumor(:,M.I.ind));
        end

    otherwise
        error("No method setup for initializing agents like %s.\n",M.setup.agent_initialization_location);

end


%% determining remaining wait times for proliferation on tumor cells (see my Onenote document in ABM>Rates--Exponential>Going backwards in time)
% old method
% prolif_rate = M.pars.prolif_rate; % assume that all mutant cells have phiD value of gammaT at start for simplicity
% u = rand(M.setup.N0,1);
% i1 = u>=1/(1+M.pars.min_prolif_wait*prolif_rate);
% i2 = ~i1;
% x = zeros(M.setup.N0,1);
% x(i1) = u(i1)*(1+M.pars.min_prolif_wait*prolif_rate)/prolif_rate - 1/prolif_rate - M.pars.min_prolif_wait;
% x(i2) = log((1+M.pars.min_prolif_wait*prolif_rate)*u(i2))/prolif_rate - M.pars.min_prolif_wait;
% M.tumor(:,M.I.proliferation_timer) = max(0,x+M.pars.min_prolif_wait); % time until next possible proliferation (days)

% new method
% r = M.pars.prolif_rate;
% w = M.pars.min_prolif_wait;
% % prob_on_clock = 1.5 + .5/(r*w) - .5*sqrt(1+6/(r*w)+1/(r*r*w*w)); % solve steady-state of x'=-(1/w)x+2ry-dx;y'=(1/w)x-ry-dy; and compute p=x/(x+y). x=# on clock, y=# off clock
% 
% I = fzero(@(I) I+1 - 2*exp(-r*I*w),0.5);
% % I = double(solve(2==(I+1)*exp(r*I*w),I));
% prob_on_clock = 1 - I; % let u = distribution on times sinces last proliferation such that \int_0^\inf u dw = # cells. then, u_t = -u_w (clock progresses) - du (apoptosis) - ruH(w-T) (proliferation after waiting T time) with u(0,t) = 2r\int_T^\inf u dw (proliferation results in 2 daughter cells with 0 wait time) and u(\inf,t)=0 (no cells have infinite time since last proliferation); then find steady-state for p = u(w,t)/N(t), N(t) = \int_0^\inf u(w,t)dw is the total number of cells at time t
% 
% on_clock_log = rand(M.setup.N0,1) < prob_on_clock;
% M.tumor(~on_clock_log,M.I.proliferation_timer) = 0;
% n_on_clock = sum(on_clock_log);
% % time_on_clock = log(1+rand(n_on_clock,1)*(exp(M.pars.apop_rate*M.pars.min_prolif_wait)-1))/M.pars.apop_rate; % assume remaining time distribution (rho) on clock satisfies PDE: rho_t = rho_tau - drho where t is the simulation time variable and tau is the remaining time; solve for steady-state: rho_t=0 to get rho(tau)=Cexp(dtau) at steady-state; normalize to a distribution and draw values using ICDF
% time_on_clock = log(1+rand(n_on_clock,1)*(exp(r*I*w)-1))/(r*I); % the solution to the PDE defined in comment to prob_on_clock = 1 - I; (above) defines a probability distribution on times since last proliferation; restricting this to [0,w] and normalize to get a probability distribution conditional on time since prolif < w and use the ICDF to get this
% M.tumor(on_clock_log,M.I.proliferation_timer) = time_on_clock;

%% set event and phase
M.tumor(:,M.I.event) = 0; % event index

%% set phase

% just_did_M = M.tumor(:,M.I.proliferation_timer)>=M.pars.min_prolif_wait-M.pars.max_dt;
% still_in_g0 = ~just_did_M & M.tumor(:,M.I.proliferation_timer)>0;
% in_g1 = ~just_did_M & ~still_in_g0;
% 
% M.tumor(just_did_M,M.I.phase) = M.val.phase_m;
% M.tumor(still_in_g0,M.I.phase) = M.val.phase_g0;
% M.tumor(in_g1,M.I.phase) = M.val.phase_g1;

p = 1./M.cycle_pars.transition_rates;
p = p/sum(p);

[M.tumor(:,M.I.phase),~] = find(diff(rand(M.setup.N0,1)<cumsum([0,p,1],2),1,2)');

M.NT = size(M.tumor,1);