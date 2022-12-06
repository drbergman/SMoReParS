function M = saveModelData(M)

time = M.t;

healthy_locations = M.save_pars.integrify(M.healthy(:,M.I.subs));
healthy_is_ha = M.healthy(:,M.I.type)==1;
healthy_is_fgfr3_mut = M.healthy(:,M.I.healthy_mut)==1;

immune_locations = M.save_pars.integrify(M.immunes(:,M.I.subs));
immune_type = int8(M.immunes(:,M.I.type));

fgfr3_concentrations = M.fgfr3.concentration;
apd1_concentrations = M.checkpoint.aPD1;

save(sprintf("data/%s/output_%08d",M.save_pars.sim_identifier,M.save_index),'-regexp','[^M]','-v7.3')


M.next_save_time = M.next_save_time + M.save_pars.dt;
M.save_index = M.save_index + 1;