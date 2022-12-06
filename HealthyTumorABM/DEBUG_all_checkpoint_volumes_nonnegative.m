function DEBUG_all_checkpoint_volumes_nonnegative(M)

assert(all(M.checkpoint.volumes_neighbors.pd1>=0,'all'))
assert(all(M.checkpoint.volumes_neighbors.pdl1>=0,'all'))