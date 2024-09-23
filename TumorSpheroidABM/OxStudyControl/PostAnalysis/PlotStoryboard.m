clearvars;

sim_id = "2403122307";
load(sprintf("../../data/sims/%s/output_constants.mat",sim_id),"grid_size")
phase_names = ["G1","S","G2","M"];
phase_size = [3,4,6,6];
phase_colors = lines(4);

f = figureOnRight;
ax = gobjects(2,2);
for day_ind = 1:4
    day = day_ind-1;
    load(sprintf("../../data/sims/%s/output_%08d.mat",sim_id,day))

    ax(day_ind) = subplot(2,2,day_ind);
    hold on


    S = gobjects(4,1);
    for i = 1:4
        phase_log = tumor_phases==i;
        S(i) = scatter(tumor_locations(phase_log,1),tumor_locations(phase_log,2),phase_size(i),phase_colors(i,:),"filled","DisplayName",phase_names(i));
    end

    axis(ax(day_ind),"square",[0.5 grid_size(1)+0.5 0.5 grid_size(2)+0.5])
    set(ax(day_ind),"XTick",[],"YTick",[]) %,"XColor","none","YColor","none")
    title(ax(day_ind),sprintf("t = %dd",day),"FontSize",8,"FontWeight","normal")
end
set(ax,"box","on")
set(ax, 'LineWidth', 0.5); % Sets border line width to 2 and color to black

ax = ax';
f.Units = "inches";

f.Position(3:4) = [2,2];

%% set margins
margin = struct("left",.01,"right",.01,"top",.1,"bottom",.01);
spacing = struct("horizontal",.01,"vertical",.09);
uniformAxisSpacing(ax,margin,spacing);

% print(f,"ABM_Snapshot_t3","-dsvg")
%%
print(f,"~/Documents/Research/SMoReParS_Comm/GlobalSensitivity/MATLAB_Figures/ABM_Storyboard","-dsvg")



