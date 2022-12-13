clearvars;

v_ha_mut = [0;0];
v_la_mut = [0.5;sqrt(3)/2];
v_la = [1;0];

v = [v_la,v_la_mut,v_ha_mut];



%%
f=figure;
scatter(v(1,:),v(2,:),'filled')
f.Position(1) = 1547;
f.Position(3:4) = [664 543];
hold on;

%%
cohort_name = "cohort_221021113614089";
path_to_cohort_file = sprintf("../data/%s/%s.mat",cohort_name,cohort_name);
load(path_to_cohort_file,"total_runs","ids")

for si = 1:total_runs
    sim_name = ids(si);

    path_to_sim_folder = sprintf("../data/%s",sim_name);

    load(sprintf("%s/output_final.mat",path_to_sim_folder),"tracked")

    a = reshape(tracked.phase_cell_hours,[],4);
    a = a(:,[1,2,4])';
    a = a./sum(a,1);


    %%
    % scatter3(v(1,:)*a,v(2,:)*a,v(3,:)*a,10,'green','filled')
    plot(v(1,:)*a,v(2,:)*a,'black',"LineWidth",1)
end