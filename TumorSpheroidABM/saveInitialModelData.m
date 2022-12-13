function M = saveInitialModelData(M)

grid_size = M.grid.size;
setup = M.setup;
pars = M.pars;
plot_pars = M.plot_pars;
flags = M.flags;
save_pars = M.save_pars;

save(sprintf("data/%s/output_constants",M.save_pars.sim_identifier),'-regexp','[^M]','-v7.3')

M = saveModelData(M);