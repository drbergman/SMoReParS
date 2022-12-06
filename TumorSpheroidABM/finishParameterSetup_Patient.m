function M = finishParameterSetup_Patient(M)

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


%% index and value stuff
M.I = buildIndices();
M.val = buildLatticeVals();

%% cell stuff
M.pars.V = M.pars.cell_width^3;
M.pars.move_rate = M.pars.move_rate_microns / M.pars.cell_width;


%% events
M = initializeEvents(M);
