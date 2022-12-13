function M = finishParameterSetup_Patient(M)

if M.setup.use_carrying_capacity_for_grid_size
    K = M.setup.carrying_capacity * M.pars.cell_width^M.setup.ndims; % carrying capacity in cubic microns
    M.setup.grid_size_microns_x = ceil(nthroot(K,M.setup.ndims));
    M.setup.grid_size_microns_y = ceil(nthroot(K/M.setup.grid_size_microns_x,M.setup.ndims-1));
    if M.setup.ndims==3
        M.setup.grid_size_microns_z = ceil(K/(M.setup.grid_size_microns_x*M.setup.grid_size_microns_x));
    end
end

if M.setup.ndims == 3
    M.setup.grid_size_microns = [M.setup.grid_size_microns_x,M.setup.grid_size_microns_y,M.setup.grid_size_microns_z];

    %% neighbor stuff
    M.pars.neighbors = allCombos(-1:1,-1:1,-1:1,'matlab');
    M.pars.neighbors(all(M.pars.neighbors==0,2),:) = []; % don't count self as neighbor
    M.pars.neighbors_VN = M.pars.neighbors;
    M.pars.neighbors_VN(sum(abs(M.pars.neighbors_VN),2)>1,:) = [];

    M.pars.neighbor_weights = 1./sqrt(sum(M.pars.neighbors.^2,2));
    M.pars.left_neighbors = M.pars.neighbors(:,1) == -1;
    M.pars.right_neighbors = M.pars.neighbors(:,1) == 1;
    M.pars.front_neighbors = M.pars.neighbors(:,2) == -1;
    M.pars.back_neighbors = M.pars.neighbors(:,2) == 1;
    M.pars.bottom_neighbors = M.pars.neighbors(:,3) == -1;
    M.pars.top_neighbors = M.pars.neighbors(:,3) == 1;

    M.pars.left_neighbors_VN_ind = find(M.pars.neighbors_VN(:,1) == -1);
    M.pars.right_neighbors_VN_ind = find(M.pars.neighbors_VN(:,1) == 1);
    M.pars.front_neighbors_VN_ind = find(M.pars.neighbors_VN(:,2) == -1);
    M.pars.back_neighbors_VN_ind = find(M.pars.neighbors_VN(:,2) == 1);
    M.pars.bottom_neighbors_VN_ind = find(M.pars.neighbors_VN(:,3) == -1);
    M.pars.top_neighbors_VN_ind = find(M.pars.neighbors_VN(:,3) == 1);

    M.pars.occmax = M.pars.occmax_3d;

else
    M.setup.grid_size_microns = [M.setup.grid_size_microns_x,M.setup.grid_size_microns_y];

    %% neighbor stuff
    M.pars.neighbors = allCombos(-1:1,-1:1,'matlab');
    M.pars.neighbors(all(M.pars.neighbors==0,2),:) = []; % don't count self as neighbor
    M.pars.neighbors_VN = M.pars.neighbors;
    M.pars.neighbors_VN(sum(abs(M.pars.neighbors_VN),2)>1,:) = [];

    M.pars.neighbor_weights = 1./sqrt(sum(M.pars.neighbors.^2,2));
    M.pars.left_neighbors = M.pars.neighbors(:,1) == -1;
    M.pars.right_neighbors = M.pars.neighbors(:,1) == 1;
    M.pars.front_neighbors = M.pars.neighbors(:,2) == -1;
    M.pars.back_neighbors = M.pars.neighbors(:,2) == 1;

    M.pars.left_neighbors_VN_ind = find(M.pars.neighbors_VN(:,1) == -1);
    M.pars.right_neighbors_VN_ind = find(M.pars.neighbors_VN(:,1) == 1);
    M.pars.front_neighbors_VN_ind = find(M.pars.neighbors_VN(:,2) == -1);
    M.pars.back_neighbors_VN_ind = find(M.pars.neighbors_VN(:,2) == 1);

    M.pars.occmax = M.pars.occmax_2d;

end

M.pars.move_rate = M.pars.move_rate_microns / M.pars.cell_width;


%% index and value stuff
M.I = buildIndices(M.setup.ndims);
M.val = buildLatticeVals();

%% events
M = initializeEvents(M);

%% issue warning if model time step and mitosis duration mismatch
if M.pars.max_dt ~= M.pars.mitosis_duration
    warning("Mismatch between (desired) dt and mitosis duration, which will result in recorded M proportions being erroneous. Post-processing can fix this, but no such script exists yet.")
end