function [out,non_border_neighbors] = getCellNeighborIndsOnGrid(cell_row,M)

out = cell_row(M.I.ind)+M.rel_pos_ind;

border_neighbors = false(8,1);

switch cell_row(M.I.subs(1))
    case 1 % then on left boundary
        border_neighbors(M.pars.left_neighbors) = true;
    case M.grid.size(1) % then on right boundary
        border_neighbors(M.pars.right_neighbors) = true;
end

switch cell_row(M.I.subs(2))
    case 1 % then on front boundary
        border_neighbors(M.pars.front_neighbors) = true;
    case M.grid.size(2) % then on back boundary
        border_neighbors(M.pars.back_neighbors) = true;
end

out(border_neighbors) = [];
non_border_neighbors = find(~border_neighbors);