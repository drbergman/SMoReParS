clearvars;

v_ha = [0;0;0];
v_ha_mut = [1;0;0];
v_la = [0.5;sqrt(3)/2;0];
v_la_mut = [.5;sqrt(3)/6;sqrt(2/3)];

v = [v_la,v_la_mut,v_ha,v_ha_mut];



%%
f=figure;
scatter3(v(1,:),v(2,:),v(3,:),'filled')
f.Position(1) = 1547;
f.Position(3:4) = [664 543];
hold on;
view(-32,9)
%%
cohort_name = "cohort_221021113614089";
path_to_cohort_file = sprintf("../data/%s/%s.mat",cohort_name,cohort_name);
load(path_to_cohort_file,"total_runs","ids")

for si = 1:total_runs
sim_name = ids(si);

path_to_sim_folder = sprintf("../data/%s",sim_name);

load(sprintf("%s/output_final.mat",path_to_sim_folder),"tracked")

a = tracked.phase_cell_days./tracked.NT;
a = reshape(a,[],4)';


%%
% scatter3(v(1,:)*a,v(2,:)*a,v(3,:)*a,10,'green','filled')
plot3(v(1,:)*a,v(2,:)*a,v(3,:)*a,'black',"LineWidth",1)
end