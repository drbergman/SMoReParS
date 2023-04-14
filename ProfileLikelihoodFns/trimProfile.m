function profile = trimProfile(profile,threshold)

[I1,I2] = getProfileInds(profile(end,:),threshold);
max_val = min(profile(end,:)) + threshold;

if I1~=1
    slope = diff(profile(:,I1+(-1:0)),1,2);
    d = (max_val - profile(end,I1)) / diff(profile(end,I1+(-1:0)));
    L = profile(:,I1) + d * slope;
else
    L = zeros(size(profile,1),0);
end
if I2~=size(profile,2)
    slope = diff(profile(:,I2+(0:1)),1,2);
    d = (max_val - profile(end,I2)) / diff(profile(end,I2+(0:1)));
    R = profile(:,I2) + d * slope;
else
    R = zeros(size(profile,1),0);
end
profile = [L,profile(:,I1:I2),R];