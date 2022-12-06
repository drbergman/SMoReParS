function M = initializeTracked(M)

M.tracked.t = 0;

M.tracked.healthy_count = M.healthy_count;
M.tracked.healthy_types = [size(M.healthy,1),0,0]; % non-cycling / cycling

M.tracked.healthy_apop = 0;
M.tracked.healthy_prolif = 0;
M.tracked.healthy_contact_inhibition = 0;

M.tracked.tumor_count = M.tumor_count;
M.tracked.tumor_types = [size(M.tumor,1),0,0]; % non-cycling / cycling

M.tracked.tumor_apop = 0;
M.tracked.tumor_prolif = 0;
M.tracked.tumor_contact_inhibition = 0;
