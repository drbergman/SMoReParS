function out = complexModel(p)

out = sum(p,1) .* (1 + 0.01 * randn(1,size(p,2)));