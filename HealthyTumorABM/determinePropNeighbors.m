function region_props = determinePropNeighbors(M)

% this function computes region_props which determines for each of three
% occupied states (nonmutant, mutant, non-healthy) in each of M.fgfr3.n_regions what
% proportion of the neighbors of those lattice sites are of any given
% (type,region) combo

% In other words, the dimensions of region_props can be understood to be
% (type of neighbor, region number of neighbor, type of center, region
% number of center) with the value there indicating the proportion of all
% neighbors to that center (type,region) combo being that particular
% (type,region) combo.

cell_loc_ind = M.healthy(:,M.I.ind);

region_props = zeros([3,M.fgfr3.n_regions,3,M.fgfr3.n_regions]);
%% get lattice indices for all neighbor sites; for those on the border, use the center site
neighbor_inds = cell_loc_ind'+M.rel_pos_ind_VN;
cells_on_left = M.healthy(:,M.I.subs(1))==1;
neighbor_inds(M.pars.left_neighbors_VN_ind,cells_on_left) = M.healthy(cells_on_left,M.I.ind);
cells_on_right = M.healthy(:,M.I.subs(1))==M.grid.size(1);
neighbor_inds(M.pars.right_neighbors_VN_ind,cells_on_right) = M.healthy(cells_on_right,M.I.ind);
cells_on_front = M.healthy(:,M.I.subs(2))==1;
neighbor_inds(M.pars.front_neighbors_VN_ind,cells_on_front) = M.healthy(cells_on_front,M.I.ind);
cells_on_back = M.healthy(:,M.I.subs(2))==M.grid.size(2);
neighbor_inds(M.pars.back_neighbors_VN_ind,cells_on_back) = M.healthy(cells_on_back,M.I.ind);
cells_on_bottom = M.healthy(:,M.I.subs(3))==1;
neighbor_inds(M.pars.bottom_neighbors_VN_ind,cells_on_bottom) = M.healthy(cells_on_bottom,M.I.ind);
cells_on_top = M.healthy(:,M.I.subs(3))==M.grid.size(3);
neighbor_inds(M.pars.top_neighbors_VN_ind,cells_on_top) = M.healthy(cells_on_top,M.I.ind);

%%
neighbor_region = M.fgfr3.regions(neighbor_inds); % region indices for all healthy neighbors
neighbor_val = full(M.L(neighbor_inds)); % get all neighbor values at once

cell_region = M.healthy(:,M.I.region_fgfr3);

myneighbors = zeros(size(region_props)); % initialize this, which is the count version of region_props
myneighbors(3,:,3,:) = M.fgfr3.neighbors_when_empty; % start by assuming all lattice sites are non-healthy, so that all non-healthy centers have all non-healthy neighbors as already computed by initializeBloodVessels; all other center (type, region) combos have no neighbors as of yet  
for center_ri = 1:M.fgfr3.n_regions % vary over the regions
    for center_ti = 1:2 % vary over the healthy types
        temp_ind = (cell_region==center_ri) & ((M.mut_log+1)==center_ti); % locations that match this region & healthy type combo
        bv_temp = neighbor_region(:,temp_ind); % region index of all the neighbors to these (type, region) combos
        myneighbors = myneighbors - accumarray([3*ones(numel(bv_temp),1),bv_temp(:),3*ones(numel(bv_temp),1),center_ri*ones(numel(bv_temp),1)],1,size(myneighbors)); % since these sites are not empty, remove these counts from this (non-healthy, center region) counts
        n_temp = neighbor_val(:,temp_ind); % healthy type of the neighbors
        for neighbor_ri = 1:M.fgfr3.n_regions
            x = bv_temp(:)==neighbor_ri; % neighbors to this (type, region) combo that are in region number (neighbor_ri)
            for neighbor_ti = 1:2
                myneighbors(neighbor_ti,neighbor_ri,center_ti,center_ri) = sum(x & n_temp(:)==M.val.healthy(neighbor_ti)); % this counts how many of the neighbors that are both in region neighbor_ri and are of type neighbor_ti (neighbors that are non-healthy counted later)
            end
            myneighbors(3,neighbor_ri,center_ti,center_ri) = sum(x) - sum(myneighbors(1:2,neighbor_ri,center_ti,center_ri)); % count all the neighbors that are in this region and subtract off those that were just found to have a healthy there to get the number that are non-healthy
            myneighbors([center_ti,3],center_ri,3,neighbor_ri) = myneighbors([center_ti,3],center_ri,3,neighbor_ri) + [1;-1]*myneighbors(3,neighbor_ri,center_ti,center_ri); % because neighbors that are non-healthy are not looped over, so update them here based on healthy lattice sites we just found; we decrement this count by the number of neighbors that had a non-healthy neighbor in this region as each of those corresponds one-to-one with non-healthy sites that have a healthy neighbor
        end
    end
end

% for each (type, region) combo, compute proportion that are of any other
% (type, region) combo
for center_ri = 1:M.fgfr3.n_regions
    for center_ti = 1:3
        if any(myneighbors(:,:,center_ti,center_ri)~=0,'all') % only do division if there is a positive number of neighbors to put in the denominator
            region_props(:,:,center_ti,center_ri) = myneighbors(:,:,center_ti,center_ri)/sum(myneighbors(:,:,center_ti,center_ri),'all');
        end
    end
end

% assert(pA.mut_prop+pA.nonmut_prop+pA.free_prop==1)