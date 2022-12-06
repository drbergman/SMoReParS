function rate_matrix = computeRateMatrix(M)

rate_matrix.healthy = zeros(M.healthy_count,2);

rate_matrix.healthy(:,1) = M.healthy_pars.prolif_rate;
rate_matrix.healthy(:,2) = M.healthy_pars.apop_rate;

rate_matrix.tumor = zeros(M.tumor_count,2);

rate_matrix.tumor(:,1) = M.tumor_pars.prolif_rate;
rate_matrix.tumor(:,2) = M.tumor_pars.apop_rate;
