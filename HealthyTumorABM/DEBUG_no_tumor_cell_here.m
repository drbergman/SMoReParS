function DEBUG_no_healthy_cell_here(M,ind)

assert(all(M.healthy(:,M.I.ind)~=ind))