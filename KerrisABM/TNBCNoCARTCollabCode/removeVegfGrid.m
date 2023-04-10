%fillVegfGrid
% This fills the vegf grid with 20 in the approximate location of the cell agents
%input:voxel grid and pos of cell 
%Version 1
%Last Revision: Kerri 3 14, 2017

function [grid] = removeVegfGrid(grid, pos)
if nargin < 2 
    error('Too few args') 
end

shiftpos = 20*pos;
startp = shiftpos-19;

grid.VEGF(startp(1):shiftpos(1),startp(2):shiftpos(2),startp(3):shiftpos(3)) = 0;
 