function M = initializeTracked(M)

M.tracked.t = 0;
M.tracked.NT = M.NT;

M.tracked.phases = sum(M.tumor(:,M.I.phase)==1:4,1);

M.tracked.tum_apop = 0;
M.tracked.tum_prolif = 0;
M.tracked.tum_contact_inhibition = 0;
M.tracked.chemo_arrest = 0;
M.tracked.arrest_recovery = 0;
