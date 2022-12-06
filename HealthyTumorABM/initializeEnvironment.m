function M = initializeEnvironment(M)

M = initializeGrid(M);

M = initializeHealthy(M);
M = initializeTumor(M);

M.L = zeros(M.V_tot,1);
M.L(M.healthy(:,M.I.ind))=M.val.healthy;


