function M = updateTracked(M)

M.tracked.t(M.i) = M.t;

M.tracked.healthy_count(M.i) = M.healthy_count;
% M.tracked.healthy_types(M.i,:) = [M.tracked.healthy_count(M.i-1),0] + M.tracked.healthy_prolif(M.i)*[-1,1];

M.tracked.tumor_count(M.i) = M.tumor_count;
% M.tracked.tumor_types(M.i,:) = [M.tracked.tumor_count(M.i-1),0] + M.tracked.tumor_prolif(M.i)*[-1,1];
