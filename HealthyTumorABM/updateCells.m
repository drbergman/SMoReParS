function [M,event_out] = updateCells(M)

[M,event_out] = eventSelection_Cells(M);

M = performEvents(M,event_out);
