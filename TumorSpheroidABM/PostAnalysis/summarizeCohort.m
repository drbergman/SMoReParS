% a script to summarize cohort data. not sure why I originally created this
% as a function...but it seems more convenient to have this as a script.
% saves mean cell counts, phase counts, ode state variable counts, and all their
% SDs.

clearvars;
% function summarizeCohort(cohort_name)
cohort_name = "cohort_2303271138";
addpath("~/Documents/MATLAB/myfunctions/")
C = load(sprintf("../data/%s/output.mat",cohort_name));

%%
for i = numel(C.ids):-1:1
    S = load(sprintf("../data/sims/%s/output_final.mat",C.ids(i)));
    phase_count(:,:,i) = S.tracked.phases;
end

count = squeeze(sum(phase_count,2));

t = S.tracked.t;
nt = length(t);

%%
count = reshape(count,[nt,size(C.ids)]);
phase_count = reshape(phase_count,[nt,4,size(C.ids)]);
ode_state_count = cat(2,sum(sliceof(phase_count,2,1:2),2),sum(sliceof(phase_count,2,3:4),2));

%%
average_count = mean(count,ndims(count));
phase_average = mean(phase_count,ndims(phase_count));
ode_state_average = mean(ode_state_count,ndims(ode_state_count));

%%
count_std = std(count,[],ndims(count));
phase_std = std(phase_count,[],ndims(phase_count));
ode_state_std = std(ode_state_count,[],ndims(ode_state_count));

%%
save(sprintf("../data/%s/summary.mat",cohort_name),"average_count","phase_average","ode_state_average","count_std","phase_std","ode_state_std","-v7.3")

%%
rmpath("~/Documents/MATLAB/myfunctions/")