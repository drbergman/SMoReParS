function [inds_first_below,inds_last_below] = getIndicesJustBelowThreshold(vals,threshold)

[inds_before_drop,inds_last_below,max_val] = findMaxCrosses(vals,threshold);
inds_first_below = inds_before_drop+1;
if vals(1)<=max_val
    inds_first_below = [1,inds_first_below];
end

if vals(end)<=max_val
    inds_last_below = [inds_last_below,length(vals)];
end
