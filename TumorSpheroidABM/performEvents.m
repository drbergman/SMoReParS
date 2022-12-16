function M = performEvents(M,in)

for ord_ind=1:length(in.active_ind)

    j = in.active_ind(ord_ind);

    switch in.events(ord_ind)
        case 1 % transition
            %%
            %             [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid_old(M.tumor(j,:),M); % neighbor indices
            phase = M.tumor(j,M.I.phase);
            if phase==M.cycle.m % then this cell proliferates
                if M.cycle_pars.dna_check(M.cycle.m) && rand() < M.cycle_pars.arrest_prob(M.cycle.m)
                    M.tumor(j,M.I.event) = 2; % mark the cell for apoptosis
                    M.tracked.chemo_arrest(M.i) = M.tracked.chemo_arrest(M.i)+1;
                    continue;
                end

                [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid(M.tumor(j,M.I.ind),M.tumor(j,M.I.subs),M); % neighbor indices
                if M.setup.ndims==3
                    if (nnz(M.L(n_ind)) + (26-length(n_ind)))<=M.pars.occmax % check how many M.pars.neighbors are occupied
                        weights = (M.L(n_ind)==0).*M.pars.neighbor_weights(non_border_neighbors);
                        %                 ind = randsample(length(n_ind),1,true,weights); % weight by whether sites are empty and by the reciprocal of the distance
                        ind = find(sum(weights)*rand() < cumsum(weights),1);
                        rel_loc = M.pars.neighbors(non_border_neighbors(ind),:); % get the relative position of the new spot
                        M.tumor(end+1,:) = 0;
                        M.tumor(end,M.I.subs) = M.tumor(j,M.I.subs)+rel_loc;  % store new array-specific locations
                        M.tumor(end,M.I.ind) = n_ind(ind); % store new array-specific locations

                        M.L(n_ind(ind)) = M.val.tum; % set value at lattice site
                        M.tumor([j,end],M.I.phase) = M.cycle.g1;

                        M.tracked.tum_prolif(M.i) = M.tracked.tum_prolif(M.i)+1;

                    else
                        M.tracked.tum_contact_inhibition(M.i) = M.tracked.tum_contact_inhibition(M.i)+1;
                        M.tumor(j,M.I.phase) = M.cycle.g1;
                    end
                else
                    if (nnz(M.L(n_ind)) + (8-length(n_ind)))<=M.pars.occmax % check how many M.pars.neighbors are occupied
                        weights = (M.L(n_ind)==0).*M.pars.neighbor_weights(non_border_neighbors);
                        %                 ind = randsample(length(n_ind),1,true,weights); % weight by whether sites are empty and by the reciprocal of the distance
                        ind = find(sum(weights)*rand() < cumsum(weights),1);
                        rel_loc = M.pars.neighbors(non_border_neighbors(ind),:); % get the relative position of the new spot
                        M.tumor(end+1,:) = 0;
                        M.tumor(end,M.I.subs) = M.tumor(j,M.I.subs)+rel_loc;  % store new array-specific locations
                        M.tumor(end,M.I.ind) = n_ind(ind); % store new array-specific locations

                        M.L(n_ind(ind)) = M.val.tum; % set value at lattice site
                        M.tumor([j,end],M.I.phase) = M.cycle.g1;

                        M.tracked.tum_prolif(M.i) = M.tracked.tum_prolif(M.i)+1;

                    else
                        M.tracked.tum_contact_inhibition(M.i) = M.tracked.tum_contact_inhibition(M.i)+1;
                        M.tumor(j,M.I.phase) = M.cycle.g1;
                    end

                end
            else % then it just moves along the transition path
                
                if M.cycle_pars.dna_check(phase) && rand() < M.cycle_pars.arrest_prob(phase)
                    M.tumor(j,M.I.event) = 2; % mark the cell for apoptosis
                    M.tracked.chemo_arrest(M.i) = M.tracked.chemo_arrest(M.i)+1;
                else
                    M.tumor(j,M.I.phase) = M.cycle.advancer(M.tumor(j,M.I.phase));
                end
            end

        case 2 % spontaneous apoptosis
            %%
            M.L(M.tumor(j,M.I.ind)) = M.val.tum_apop;
            M.tracked.tum_apop(M.i) = M.tracked.tum_apop(M.i)+1;

        case 3 % random movement
            %%
            [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid(M.tumor(j,M.I.ind),M.tumor(j,M.I.subs),M); % neighbor indices
            if any(M.L(n_ind)==0)
                weights = (M.L(n_ind)==0).*M.pars.neighbor_weights(non_border_neighbors);
                ind = find(sum(weights)*rand() < cumsum(weights),1);
                rel_loc = M.pars.neighbors(non_border_neighbors(ind),:); % get the relative position of the new spot
                M.tumor(j,M.I.subs) = M.tumor(j,M.I.subs)+rel_loc;  % store new array-specific locations

                M.L(M.tumor(j,M.I.ind)) = 0;

                M.tumor(j,M.I.ind) = n_ind(ind); % store new array-specific locations
                M.L(n_ind(ind)) = M.val.tum; % set value at lattice site based on original cell (this should work instead of above commented out line)
            end

        case 4 % spontaneous apoptosis
            %%
            error("no longer have this")
            M.L(M.tumor(j,M.I.ind)) = M.val.tum_apop;
            M.tracked.chemo_arrest(M.i) = M.tracked.chemo_arrest(M.i)+1;
            M.tumor(j,M.I.event) = 2; % mark it as apoptotic to be removed later

        otherwise
            %%
            error('should not do nothing')

    end % end of switch
end % end of j for