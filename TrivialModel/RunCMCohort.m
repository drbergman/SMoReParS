% Run a cohort of the complex model.

clearvars

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions

a = 1:3;
b = 4:7;

nsamps = 100;

[P(:,:,1),P(:,:,2)] = ndgrid(a,b);

P = reshape(P,[],2)';

for si = 1:nsamps
out(:,:,si) = complexModel(P);
end
out = reshape(out,length(a),length(b),nsamps);

out_avg = mean(out,3);
out_std = std(out,[],3);
for i = 1:length(a)
    for j = 1:length(b)
        D(1,i,j).A = out_avg(i,j);
        D(1,i,j).S = out_std(i,j);
    end
end

C = {[]};
t = 0;

par_names = ["a","b"];

vals = {a',b'};

cohort_folder_spec = "data/cohort_%d";
cohort_num = next_version_number(cohort_folder_spec);
cohort_folder_string = sprintf(cohort_folder_spec,cohort_num);

cohort_size = [length(a),length(b)];
mkdir(cohort_folder_string)
save(sprintf("%s/summary.mat",cohort_folder_string),"cohort_size","nsamps","D","par_names","vals","C","t")