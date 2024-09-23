clearvars;

sim_id = "2403122307";
day = 2;
load(sprintf("../../data/sims/%s/output_constants.mat",sim_id),"grid_size")
load(sprintf("../../data/sims/%s/output_%08d.mat",sim_id,day))

f = figure;
ax = gca;
hold on

phase_names = ["G1","S","G2","M"];
phase_size = [3,4,6,6];
phase_colors = lines(4);

S = gobjects(4,1);
for i = 1:4
    phase_log = tumor_phases==i;
    S(i) = scatter(tumor_locations(phase_log,1),tumor_locations(phase_log,2),phase_size(i),phase_colors(i,:),"filled","DisplayName",phase_names(i));
end

axis(ax,"square",[0.5 grid_size(1)+0.5 0.5 grid_size(2)+0.5])
set(ax,"XTick",[],"YTick",[],"XColor","none","YColor","none")

L = legend(ax,"Location","bestoutside","FontSize",6);

f.Units = "inches";
f.Position(3:4) = [2,1.25];
ax.Units = "inches";
ax.Position = [0.15 0.15 1 1];

% print(f,"ABM_Snapshot_t3","-dsvg")
print(f,sprintf("~/Documents/Research/SMoReParS_Comm/GlobalSensitivity/MATLAB_Figures/ABM_Snapshot_t%d",day),"-dsvg")



