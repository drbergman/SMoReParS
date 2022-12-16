function rate_matrix = computeRateMatrix(M)

rate_matrix = zeros(M.NT,3);

rate_matrix(:,1) = M.cycle_pars.transition_rates(M.tumor(:,M.I.phase));
rate_matrix(:,2) = M.pars.apop_rate;
rate_matrix(:,3) = M.pars.move_rate;
