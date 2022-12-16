function M = updateTracked(M)

M.tracked.t(M.i) = M.t;

M.tracked.NT(M.i) = M.NT;

M.tracked.phases(M.i,:) = sum(M.tumor(:,M.I.phase)==1:4,1);
