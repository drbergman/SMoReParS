% This script will convert a profile of LMS using d and Delta d into one
% using the Hill parameters (delta and the ec50)

clearvars;

addpath("../ODEFitting/")

file_name = "Profiles_SMFromData_LMS_bounded_clean";

load("data/" + file_name,"profiles");

I = [8,9]; 
b = 4;
% if ~input("Are you sure that they are numbers 8 and 9?")
%     error("Make sure this is correct.")
% end
for i = 1:numel(profiles)
    profiles{i} = updateProfile(profiles{i},I,b);
end

save("data/" + file_name + "_converted","profiles");


%% reset path
rmpath("../ODEFitting/")


function profile = updateProfile(profile,I,b)

d = profile(I(1),:);
delta_d = profile(I(2),:);

[profile(I(1),:),profile(I(2),:)] = computeHillParameters(d,delta_d,b);

end
