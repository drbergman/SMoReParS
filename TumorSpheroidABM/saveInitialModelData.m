function M = saveInitialModelData(M)

grid_size = M.grid.size;
setup = M.setup;
pars = M.pars;
cycle_pars = M.cycle_pars;
plot_pars = M.plot_pars;
flags = M.flags;
save_pars = M.save_pars;

save(sprintf("data/sims/%s/output_constants",M.save_pars.sim_identifier),'-regexp','[^M]','-v7.3')

if M.save_pars.dt < Inf
    M = saveModelData(M);
end
