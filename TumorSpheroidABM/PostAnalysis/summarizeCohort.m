function summarizeCohort(cohort_name)

addpath("~/Documents/MATLAB/myfunctions/")
C = load(sprintf("../data/%s/output.mat",cohort_name));

%%
for i = numel(C.ids):-1:1
    S = load(sprintf("../data/sims/%s/output_final.mat",C.ids(i)));
    count(:,i) = S.tracked.NT;
    phase_count(:,:,i) = S.tracked.phases;
end

t = S.tracked.t;
nt = length(t);

%%
count = reshape(count,[nt,size(C.ids)]);
phase_count = reshape(phase_count,[nt,4,size(C.ids)]);

%%
average_count = mean(count,ndims(count));
phase_average = mean(phase_count,ndims(phase_count));

%%
count_std = std(count,[],ndims(count));
phase_std = std(phase_count,[],ndims(phase_count));

%%
save(sprintf("../data/%s/summary.mat",cohort_name),"average_count","phase_average","count_std","phase_std","-v7.3")

%%
rmpath("~/Documents/MATLAB/myfunctions/")