% a script to collate all ABM sims from a single cohort.

clearvars;
cohort_name = "cohort_2306062212";
addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
load(sprintf("../../data/%s/output.mat",cohort_name),"ids","nsamps_per_condition","cohort_size","lattice_parameters");
nsamps_per_parameter_vector = nsamps_per_condition;
n_conditions = 1;
C = {[]};
vals = {lattice_parameters.values};
%%
t = [0;10;24;36;48;72]/24;
nt = length(t);

for i = numel(ids):-1:1
    S = load(sprintf("../../data/sims/%s/output_final.mat",ids(i)));
    t_abm = S.tracked.t;
    t_abm = round(1440*t_abm)/1440; % to make sure that the last time point is actually 3 days (not 3-eps() days)
    temp = S.tracked.phases;
    temp = interp1(t_abm,temp,t);
    phase_count(:,:,i) = temp;
    n_phases(i) = size(S.tracked.phases,2);
end

n_phases = unique(n_phases);
if length(n_phases)>1
    error("different numbers of phases in these sims")
end

phase_count = reshape(phase_count,[nt,n_phases,cohort_size,nsamps_per_condition]);
phase_count_dims = ["time","phase (G1,S,G2,M,G1 arrest,S arrest,G2 arrest, M arrest)","dose","varied parameters (x5)","samples"];

save(sprintf("%s_collated.mat",cohort_name),"t","phase_count","lattice_parameters","phase_count_dims","-v7.3")
