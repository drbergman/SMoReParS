clearvars;

addpath("../..")
addpath("~/Documents/MATLAB/myfunctions/")
cohort_id = "cohort_230124175743017";

load(sprintf("../../data/%s/output.mat",cohort_id),"ids","lattice_parameters","nsamps_per_condition")

ncohorts = numel(ids)/nsamps_per_condition;
npars = numel(lattice_parameters);

cycle = buildCycle();

ids = reshape(ids,[],nsamps_per_condition);
I = randi(size(ids,1));
f=figure;
ax = gca;
hold on;
t = zeros(289,nsamps_per_condition);
N = zeros(289,nsamps_per_condition,4); % may need to change this once the arrested compartment is finalized
for si = 1:nsamps_per_condition
    T(si) = load(sprintf("../../data/sims/%s/output_final.mat",ids(I,si)),"tracked");
    t(:,si) = T(si).tracked.t;
    N(:,si,:) = T(si).tracked.phases; 
end

colors = lines(4);
names = ["G_1","S","G_2","M"];
l = gobjects(4,1);
opts.min_val = 0;
for pi = 1:4 % phase index
    opts.Color = colors(pi,:);
    opts.DisplayName = names(pi);
    [~,l(pi)] = patchPlot(ax,t,N(:,:,pi),opts);
end

legend(ax,l,"location","best")
xlabel("Time (d)")
ylabel("Count")
set(ax,"FontSize",16)

file_name = "TimeSeriesExample";
savefig(f,sprintf("figures/fig/%s",file_name))
print(f,sprintf("figures/png/%s",file_name),"-dpng")

rmpath("../..")
