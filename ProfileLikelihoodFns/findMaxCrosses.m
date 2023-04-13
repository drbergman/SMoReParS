function [I1,I2,max_val] = findMaxCrosses(vals,threshold)

min_val = min(vals);
max_val = min_val+threshold;
changes = diff(vals>max_val);

I1 = find(changes==-1); % times the profile dipped below the threshold, i.e. where the min should par value should be
I2 = find(changes==1); % times the profile rose above the threshold, i.e. where the max should par value should be
