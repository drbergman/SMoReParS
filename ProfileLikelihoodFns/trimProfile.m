function profile = trimProfile(profile,threshold)

[par_min,par_max] = getProfileBounds(profile,threshold);
max_val = min(profile(end,:)) + threshold;
I = profile(1,:) >= par_min & profile(1,:) <= par_max;
profile = profile(:,I);
if profile(1,1)>par_min
    profile = [[par_min;max_val],profile];
end
if profile(1,end)<par_max
    profile = [profile,[par_max;max_val]];
end