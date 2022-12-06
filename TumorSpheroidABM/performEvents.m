function M = performEvents(M,in)

for ord_ind=1:length(in.active_ind)

    j = in.active_ind(ord_ind);

    switch in.events(ord_ind)
        case 1 % proliferation
            %%
            %             [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid_old(M.tumor(j,:),M); % neighbor indices
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

                    M.tracked.tumor_types(M.i,[2,3]) = M.tracked.tumor_types(M.i,[2,3]) + (min(M.pars.mitosis_duration,M.dt-M.tumor(j,M.I.proliferation_timer))/M.dt) * [-1,1];

                    M.tumor([j,end],M.I.proliferation_timer) = M.pars.min_prolif_wait + in.time_to_event(ord_ind); % must wait M.pars.min_prolif_wait days before proliferating; assume the proliferation happened sometime in the interval [M.tumor(j,M.I.proliferation_timer),M.dt] so some progress towards next prolif has happened (will subtract off M.dt with all cells in simForward)

                    M.tracked.tum_prolif(M.i) = M.tracked.tum_prolif(M.i)+1;

                else
                    M.tumor(j,M.I.proliferation_timer) = 0; % if not enough empty space, then allow this cell to try proliferating again
                    M.tracked.tum_contact_inhibition(M.i) = M.tracked.tum_contact_inhibition(M.i)+1;
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

                    M.tracked.tumor_types(M.i,[2,3]) = M.tracked.tumor_types(M.i,[2,3]) + (min(M.pars.mitosis_duration,M.dt-M.tumor(j,M.I.proliferation_timer))/M.dt) * [-1,1];

                    M.tumor([j,end],M.I.proliferation_timer) = M.pars.min_prolif_wait + in.time_to_event(ord_ind); % must wait M.pars.min_prolif_wait days before proliferating; assume the proliferation happened sometime in the interval [M.tumor(j,M.I.proliferation_timer),M.dt] so some progress towards next prolif has happened (will subtract off M.dt with all cells in simForward)

                    M.tracked.tum_prolif(M.i) = M.tracked.tum_prolif(M.i)+1;

                else
                    M.tumor(j,M.I.proliferation_timer) = 0; % if not enough empty space, then allow this cell to try proliferating again
                    M.tracked.tum_contact_inhibition(M.i) = M.tracked.tum_contact_inhibition(M.i)+1;
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

        otherwise
            %%
            error('should not do nothing')

    end % end of switch
end % end of j for