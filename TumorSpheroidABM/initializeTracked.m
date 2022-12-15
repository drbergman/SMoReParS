function M = initializeTracked(M)

M.tracked.t = 0;
M.tracked.NT = M.NT;

% just_did_M = M.tumor(:,M.I.proliferation_timer)>=M.pars.min_prolif_wait-M.pars.max_dt;
% still_in_g0 = ~just_did_M & M.tumor(:,M.I.proliferation_timer)>0;
% in_g1 = ~just_did_M & ~still_in_g0;
% 
% phase_cell_days = [sum(still_in_g0),sum(in_g1),sum(just_did_M)] * M.pars.max_dt; % G1, G2, M

phase_cell_days = [0,0,0];
phase_cell_days(M.val.phase_g0) = M.pars.max_dt * sum(M.tumor(:,M.I.phase)==M.val.phase_g0,1);
phase_cell_days(M.val.phase_g1) = M.pars.max_dt * sum(M.tumor(:,M.I.phase)==M.val.phase_g1,1);
phase_cell_days(M.val.phase_m) = M.pars.max_dt * sum(M.tumor(:,M.I.phase)==M.val.phase_m,1) * 0.5; % multiply be 1/2 to account for half of these cells not being present in "previous" time step because they are new daughter cells, i.e. protect against double counting this time which makes it agree with how I count things later

M.tracked.phase_cell_days = phase_cell_days;
M.tracked.simple_types = [sum(phase_cell_days([M.val.phase_g0,M.val.phase_g1])),phase_cell_days(M.val.phase_m)]; % did not just proliferate, just proliferated

M.tracked.tum_apop = 0;
M.tracked.tum_prolif = 0;
M.tracked.tum_contact_inhibition = 0;
M.tracked.tum_chemo_death = 0;
