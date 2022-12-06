function M = recruitImmune(M)

if M.flags.fgfr3_affects_immune_recruit
    rate = M.dt*M.immune_pars.immune_recruit_rate*M.healthy_count*(M.immune_pars.min_imm_recruit_prop+(1-M.immune_pars.min_imm_recruit_prop)/(1+mean(M.tracked.phiD_mean(M.i))/M.immune_pars.min_imm_recruit_prop_ec50));
else
    rate = M.dt*M.immune_pars.immune_recruit_rate*M.healthy_count;
end
n_newI = poissrnd(rate);
if n_newI>0

    M = placeImmune(M,n_newI);

%     new_immune_neighbors_inds = immunes(1:n_newI,ind_ind)'+rel_pos_ind_VN;
%     healthy_neighbors_log = reshape(any(reshape(L(new_immune_neighbors_inds),[],1)==tum_vals,2),6,[]); % I probably need to do some reshaping and stuff to get this check right

end
