function profile = trimProfile(profile,threshold)

[inds_before_drop,inds_before_rise] = getProfileBoundaryInds(profile(end,:),threshold);
max_val = min(profile(end,:)) + threshold;

if inds_before_drop~=1
    slope = diff(profile(:,inds_before_drop+(-1:0)),1,2);
    d = (max_val - profile(end,inds_before_drop)) / diff(profile(end,inds_before_drop+(-1:0)));
    L = profile(:,inds_before_drop) + d * slope;
else
    L = zeros(size(profile,1),0);
end
if inds_before_rise~=size(profile,2)
    slope = diff(profile(:,inds_before_rise+(0:1)),1,2);
    d = (max_val - profile(end,inds_before_rise)) / diff(profile(end,inds_before_rise+(0:1)));
    R = profile(:,inds_before_rise) + d * slope;
else
    R = zeros(size(profile,1),0);
end
profile = [L,profile(:,inds_before_drop:inds_before_rise),R];