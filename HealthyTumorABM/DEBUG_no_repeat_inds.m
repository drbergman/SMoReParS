function DEBUG_no_repeat_inds(M)

assert(length(unique(M.healthy(:,M.I.ind)))==size(M.healthy,1))