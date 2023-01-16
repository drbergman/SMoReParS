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

%% set event and phase
M.tumor(:,M.I.event) = 0; % event index

%% set phase
p = 1./M.cycle_pars.transition_rates;
p = p/sum(p);

[M.tumor(:,M.I.phase),~] = find(diff(rand(M.setup.N0,1)<cumsum([0,p,1],2),1,2)');

M.NT = size(M.tumor,1);