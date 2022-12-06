function M = performEvents(M,in)

for ord_ind=1:length(in.active_ind)

    j = in.active_ind(ord_ind);
    cell_type = in.cell_type(ord_ind);
    cell_pars_field = sprintf("%s_pars",cell_type);

    switch in.events(ord_ind)
        case 1 % proliferation
            %%
            [n_ind,non_border_neighbors] = getCellNeighborIndsOnGrid(M.(cell_type)(j,:),M); % neighbor indices
            if (nnz(M.L(n_ind)) + (8-length(n_ind)))<=M.(cell_pars_field).occmax % check how many M.pars.neighbors are occupied
                ind = randsample(length(n_ind),1,true,(M.L(n_ind)==0)./M.pars.neighbor_weights(non_border_neighbors)); % weight by whether sites are empty and by the reciprocal of the distance
                rel_loc = M.pars.neighbors(non_border_neighbors(ind),:); % get the relative position of the new spot
                M.(cell_type)(end+1,:) = 0;
                M.(cell_type)(end,M.I.subs) = M.(cell_type)(j,M.I.subs)+rel_loc;  % store new array-specific locations
                M.(cell_type)(end,M.I.ind) = n_ind(ind); % store new array-specific locations

                M.L(n_ind(ind)) = M.L(M.(cell_type)(j,M.I.ind)); % set value at lattice site based on original cell (this should work instead of above commented out line)

                type_name = sprintf("%s_types",cell_type);
                M.tracked.(type_name)(M.i,[2,3]) = M.tracked.(type_name)(M.i,[2,3]) + (min(M.pars.mitosis_duration,M.dt-M.(cell_type)(j,M.I.proliferation_timer))/M.dt) * [-1,1];
                M.(cell_type)([j,end],M.I.proliferation_timer) = M.(cell_pars_field).min_prolif_wait + in.time_to_event(ord_ind); % must wait M.(cell_pars_field).min_prolif_wait days before proliferating; assume the proliferation happened sometime in the interval [M.(cell_type)(j,M.I.proliferation_timer),M.dt] so some progress towards next prolif has happened (will subtract off M.dt with all cells in simForward)
    
                prolf_name = sprintf("%s_prolif",cell_type);
                M.tracked.(prolf_name)(M.i) = M.tracked.(prolf_name)(M.i)+1;
            else
                M.(cell_type)(j,M.I.proliferation_timer) = 0; % if not enough empty space, then allow this cell to try proliferating again
                
                contact_name = sprintf("%s_contact_inhibition",cell_type);
                M.tracked.(contact_name)(M.i) = M.tracked.(contact_name)(M.i)+1;
            end
            
        case 2 % spontaneous apoptosis
            %%
            apop_name = sprintf("%s_apop",cell_type);
            M.L(M.(cell_type)(j,M.I.ind)) = M.val.(apop_name);
            M.tracked.(apop_name)(M.i) = M.tracked.(apop_name)(M.i)+1;
        otherwise
            %%
            error('should not do nothing')
            
    end % end of switch
end % end of j for