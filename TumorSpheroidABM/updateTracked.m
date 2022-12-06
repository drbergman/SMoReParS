function M = updateTracked(M)

M.tracked.t(M.i) = M.t;

M.tracked.NT(M.i) = M.NT;

M.tracked.simple_types(M.i,:) = [M.tracked.NT(M.i-1),0] + M.tracked.tum_prolif(M.i) * [-1,1];
