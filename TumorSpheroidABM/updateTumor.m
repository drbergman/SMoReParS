function [M,event_out] = updateTumor(M)

[M,event_out] = eventSelection(M);

M = performEvents(M,event_out);
