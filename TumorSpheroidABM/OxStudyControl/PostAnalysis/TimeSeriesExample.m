clearvars;

addpath("../..")
addpath("~/Documents/MATLAB/myfunctions/")
cohort_id = "cohort_2403130846";

load(sprintf("../../data/%s/output.mat",cohort_id),"ids","lattice_parameters","nsamps_per_condition")

ncohorts = numel(ids)/nsamps_per_condition;
npars = numel(lattice_parameters);

cycle = buildCycle();

ids = reshape(ids,[],nsamps_per_condition);
I = randi(size(ids,1));
f=figure;
ax = gca;
hold on;
t = zeros(73,nsamps_per_condition);
N = zeros(73,nsamps_per_condition,8); % may need to change this once the arrested compartment is finalized
for si = 1:nsamps_per_condition
    T(si) = load(sprintf("../../data/sims/%s/output_final.mat",ids(I,si)),"tracked");
    t(:,si) = T(si).tracked.t;
    N(:,si,:) = T(si).tracked.phases; 
end
N = N(:,:,1:4); % these are control sims, so no arrested compartments

colors = lines(4);
names = ["G_1","S","G_2","M"];
l = gobjects(4,1);
opts.min_val = 0;
for pi = 1:4 % phase index
    [~,l(pi)] = patchPlot(ax,t,N(:,:,pi),Color=colors(pi,:),DisplayName=names(pi),min_val=0,LineWidth=0.5);
end
l = flip(l);
[~,l(end+1)] = patchPlot(ax,t,sum(N,3),Color=[0 0 0],DisplayName="Total",min_val=0,LineWidth=0.5);
yline(500,"LineStyle","--","LineWidth",0.5) % the carrying capacity for these sims

% legend(ax,flip(l),"location","best")
xlabel("Time (d)")
ylabel("Count")

%% set up figure
f.Units = "inches";
f.Position(3) = 2;
f.Position(4) = 0.75;
set(ax,"FontSize",8)

%%

file_name = "TimeSeriesExample";
% savefig(f,sprintf("figures/fig/%s",file_name))
% print(f,sprintf("figures/png/%s",file_name),"-dpng")
print('-r300',f,sprintf("~/Documents/Research/SMoReParS_Comm/GlobalSensitivity/MATLAB_Figures/%s",file_name),"-dpng")

%% combine to match sm
N_combined = reshape(N,73,nsamps_per_condition,2,2);
N_combined = sum(N_combined,3);
N_combined = squeeze(N_combined);

colors = lines(4);
colors = colors([1,3],:);
names = ["G_1/S","G_2/M"];
l = gobjects(2,1);
opts.min_val = 0;

f=figure;
ax = gca;
hold on;

for pi = 1:2 % phase index
    [~,l(pi)] = patchPlot(ax,t,N_combined(:,:,pi),Color=colors(pi,:),DisplayName=names(pi),min_val=0,LineWidth=0.5);
end
[~,l(end+1)] = patchPlot(ax,t,sum(N_combined,3),Color=[0 0 0],DisplayName="Total",min_val=0,LineWidth=0.5);
yline(500,"LineStyle","--","LineWidth",0.5) % the carrying capacity for these sims

% legend(ax,l,"location","best")
xlabel("Time (d)")
ylabel("Count")

%% set up figure
f.Units = "inches";
f.Position(3) = 1;
f.Position(4) = 0.75;
set(ax,"FontSize",8)

%%

file_name = "TimeSeriesExample_Combined";
% savefig(f,sprintf("figures/fig/%s",file_name))
% print(f,sprintf("figures/png/%s",file_name),"-dpng")
print('-r300',f,sprintf("~/Documents/Research/SMoReParS_Comm/GlobalSensitivity/MATLAB_Figures/%s",file_name),"-dpng")


rmpath("../..")
