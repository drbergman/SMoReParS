function rate_matrix = computeRateMatrix(M)

rate_matrix = zeros(M.NT,2);

rate_matrix(:,1) = M.pars.prolif_rate;
rate_matrix(:,2) = M.pars.apop_rate;
rate_matrix(:,3) = M.pars.move_rate;
