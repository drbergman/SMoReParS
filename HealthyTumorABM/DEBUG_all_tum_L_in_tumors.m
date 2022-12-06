function DEBUG_all_tum_L_in_healthys(M)

L_ind = find(any(M.L==M.val.healthy,2));
T_ind = sort(M.healthy(:,M.I.ind));
assert(isempty(setdiff(L_ind,T_ind)))