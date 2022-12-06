function DEBUG_all_healthys_on_L(M)

assert(all(M.L(M.healthy(:,M.I.ind))~=0))


assert(all(any(M.L(M.healthy(:,M.I.ind))==M.val.healthy,2)))