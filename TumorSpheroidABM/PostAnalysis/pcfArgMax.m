function out = pcfArgMax(pcf)

% computes the radius r at which the highest value of g occurs for each
% time
% Input variables:
%     pcf: nr x nt array where rows correspond to radii and columns
%     correspond to times

[~,out] = max(pcf,[],1);

