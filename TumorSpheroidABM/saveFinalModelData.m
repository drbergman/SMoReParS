function M = saveFinalModelData(M)

if M.save_pars.save_model_state
    load(sprintf("data/%s/output_%08d.mat",M.save_pars.sim_identifier,M.save_index-1),"time");
    if M.t>time % make new save
        M = saveModelData(M);
    end
end

tracked = M.tracked;
save(sprintf("data/%s/output_final",M.save_pars.sim_identifier),"tracked",'-v7.3')