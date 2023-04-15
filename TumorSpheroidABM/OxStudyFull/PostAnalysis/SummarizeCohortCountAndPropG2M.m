% a script to summarize cohort data just with count and proportion in g2/m.

clearvars;
% function summarizeCohort(cohort_name)
cohort_name = "cohort_2303301105";
addpath("~/Documents/MATLAB/myfunctions/")
C = load(sprintf("../data/%s/output.mat",cohort_name));

%%
S = load(sprintf("../data/sims/%s/output_final.mat",C.ids(1)));
phase_count = zeros(numel(S.tracked.t),4,numel(C.ids));
for i = 1:numel(C.ids)
    S = load(sprintf("../data/sims/%s/output_final.mat",C.ids(i)));
    phase_count(:,:,i) = S.tracked.phases;
end

S = load(sprintf("../data/sims/%s/output_final.mat",C.ids(1)));
count = squeeze(sum(phase_count,2));

t = S.tracked.t;
nt = length(t);

%%
count = reshape(count,[nt,size(C.ids)]);
count_average = mean(count,ndims(count));
count_std = std(count,[],ndims(count));

phase_count = reshape(phase_count,[nt,2,2,size(C.ids)]); % [time, [phase], sample] where [phase] has the 4 phases in this format [G1,G2;S,M] so summing along columns gives the two halves G1/S and G2/M
ode_state_count = squeeze(sum(phase_count,2));
state2_prop = squeeze(sliceof(ode_state_count,2,2))./count;
state2_prop_mean = mean(state2_prop,ndims(state2_prop));
state2_prop_std = std(state2_prop,[],ndims(state2_prop));

%%
filename = sprintf("../data/%s/summary.mat",cohort_name);
if exist(filename,"file")
    save(filename,"count_*","state2_prop_*","-append")
else
    save(filename,"count_*","state2_prop_*","-v7.3")
end
