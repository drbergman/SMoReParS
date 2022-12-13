function M = initializeTracked(M)

M.tracked.t = 0;
M.tracked.NT = M.NT;

just_did_M = M.tumor(:,M.I.proliferation_timer)>=M.pars.min_prolif_wait-M.pars.max_dt;
still_in_g0 = ~just_did_M & M.tumor(:,M.I.proliferation_timer)>0;
in_g1 = ~just_did_M & ~still_in_g0;

phase_cell_hours = [sum(still_in_g0),sum(in_g1),sum(just_did_M)] * M.pars.max_dt; % G1, G2, M

M.tracked.phase_cell_hours = phase_cell_hours;
M.tracked.simple_types = [sum(phase_cell_hours(1:2)),phase_cell_hours(3)]; % did not just proliferate, just proliferated

M.tracked.tum_apop = 0;
M.tracked.tum_prolif = 0;
M.tracked.tum_contact_inhibition = 0;
M.tracked.tum_chemo_death = 0;
