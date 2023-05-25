function rate_matrix = computeRateMatrix(M)

rate_matrix = zeros(M.NT,3);

rate_matrix(:,1) = M.cycle_pars.transition_rates(M.tumor(:,M.I.phase));
if M.flags.apop_only_for_arrested
    rate_matrix(M.tumor(:,M.I.is_arrested)==true,2) = M.pars.apop_rate + M.chemo_pars.apop_rate;
else
    rate_matrix(:,2) = M.pars.apop_rate + (M.tumor(:,M.I.is_arrested)==true)*M.chemo_pars.apop_rate;
end
rate_matrix(~M.tumor(:,M.I.is_arrested),3) = M.pars.move_rate;
