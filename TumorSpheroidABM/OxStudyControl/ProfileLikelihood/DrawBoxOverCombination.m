clearvars;

load("data/Profiles_SMFromABM_New_clean.mat")

save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig";"png"];
save_fig_opts.resolution = '-r1200';

n_sm_pars = size(profiles,1);
profiles = reshape(profiles,n_sm_pars,[]);
n_abm_pars = size(profiles,2);

%% accept box
face_color = [.1 .1 .1];
f = openfig("figures/fig/ParTriplesAll_SMFitToData_New.fig");
f.Name = "AcceptABMParameterBox";
ax = gca;
set(ax.Children,"LineWidth",0.5,"Marker","none","LineStyle","--")
view(22,42)
% legend("AutoUpdate","off")
ax.Children = flip(ax.Children);
for i = 1:3
    xx{i} = ax.Children(i).XData;
    yy{i} = ax.Children(i).YData;
    zz{i} = ax.Children(i).ZData;
    colors(i,:) = ax.Children(i).Color;
    disp_names{i} = ax.Children(i).DisplayName;
end

delete(ax.Children)

hold on
abm_par_ind = 2;
bounds = zeros(2,n_sm_pars);
for sm_par_ind = 1:n_sm_pars
    bounds(:,sm_par_ind) = profiles{sm_par_ind,abm_par_ind}(sm_par_ind,[1,end]);
end
% faces with constant z
p(1) = patch([bounds(:,1);flip(bounds(:,1))],repelem(bounds(:,2),2,1),bounds(1,3)*ones(4,1),face_color,"FaceAlpha",0.25);
p(2) = patch([bounds(:,1);flip(bounds(:,1))],repelem(bounds(:,2),2,1),bounds(2,3)*ones(4,1),face_color,"FaceAlpha",0.25);
plot3([0,0],[0,0],bounds(:,3),"LineWidth",1,"Color","k")

% faces with constant y
p(3) = patch([bounds(:,1);flip(bounds(:,1))],bounds(1,2)*ones(4,1),repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
p(4) = patch([bounds(:,1);flip(bounds(:,1))],bounds(2,2)*ones(4,1),repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
plot3(3*ones(1,2),bounds(:,2),[0,0],"LineWidth",1,"Color","k")

% faces with constant x
p(5) = patch(bounds(1,1)*ones(4,1),[bounds(:,2);flip(bounds(:,2))],repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
p(6) = patch(bounds(2,1)*ones(4,1),[bounds(:,2);flip(bounds(:,2))],repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
plot3(bounds(:,1),[0,0],[0,0],"LineWidth",1,"Color","k")


for i = 1:3
    P1 = interp1(yy{i},[xx{i}',yy{i}',zz{i}'],bounds(2,2));
    I1 = yy{i}>=bounds(2,2);
    c{1} = [xx{i}(I1),P1(1)];
    c{2} = [yy{i}(I1),P1(2)];
    c{3} = [zz{i}(I1),P1(3)];
    [~,order] = sort(c{i});
    plot3(c{1}(order),c{2}(order),c{3}(order),"LineStyle","--","LineWidth",0.5,"Color",colors(i,:));
    
    P2 = interp1(xx{i},[xx{i}',yy{i}',zz{i}'],bounds(2,1));
    I3 = xx{i}>=bounds(2,1);
    c{1} = [xx{i}(I3),P2(1)];
    c{2} = [yy{i}(I3),P2(2)];
    c{3} = [zz{i}(I3),P2(3)];
    [~,order] = sort(c{i});
    plot3(c{1}(order),c{2}(order),c{3}(order),"LineStyle","--","LineWidth",0.5,"Color",colors(i,:));
    
    I2 = ~I1 & ~I3;
    c{1} = [xx{i}(I2),P1(1),P2(1)];
    c{2} = [yy{i}(I2),P1(2),P2(2)];
    c{3} = [zz{i}(I2),P1(3),P2(3)];
    [~,order] = sort(c{i});
    ll(i) = plot3(c{1}(order),c{2}(order),c{3}(order),"LineStyle","-","LineWidth",0.5,"Color",colors(i,:),"DisplayName",disp_names{i});
end

xlim([0 3])
ylim([0 10])
zlim([0 2000])
% legend(ax,ll,"FontSize",20)
legend(ax,"off")
grid on
ax.FontSize = 8;

saveFigures(f,save_fig_opts)

%% reject box
face_color = [.1 .1 .1];
f = openfig("figures/fig/ParTriplesAll_SMFitToData_New.fig");
f.Name = "RejectABMParameterBox";
ax = gca;
set(ax.Children,"LineWidth",0.5,"Marker","none","LineStyle","--")
view(22,42)
% legend("AutoUpdate","off")
ax.Children = flip(ax.Children);
for i = 1:3
    xx{i} = ax.Children(i).XData;
    yy{i} = ax.Children(i).YData;
    zz{i} = ax.Children(i).ZData;
    colors(i,:) = ax.Children(i).Color;
    disp_names{i} = ax.Children(i).DisplayName;
end

hold on

ll = gobjects(3,1);
for i = 1:3
    ll(i) = plot3(0,0,0,"Color",colors(i,:),"LineWidth",0.5,"DisplayName",disp_names{i});
end

abm_par_ind = 12;
bounds = zeros(2,n_sm_pars);
for sm_par_ind = 1:n_sm_pars
    bounds(:,sm_par_ind) = profiles{sm_par_ind,abm_par_ind}(sm_par_ind,[1,end]);
end
% faces with constant z
p(1) = patch([bounds(:,1);flip(bounds(:,1))],repelem(bounds(:,2),2,1),bounds(1,3)*ones(4,1),face_color,"FaceAlpha",0.25);
p(2) = patch([bounds(:,1);flip(bounds(:,1))],repelem(bounds(:,2),2,1),bounds(2,3)*ones(4,1),face_color,"FaceAlpha",0.25);
plot3([0,0],[0,0],bounds(:,3),"LineWidth",1,"Color","k")

% faces with constant y
p(3) = patch([bounds(:,1);flip(bounds(:,1))],bounds(1,2)*ones(4,1),repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
p(4) = patch([bounds(:,1);flip(bounds(:,1))],bounds(2,2)*ones(4,1),repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
plot3(3*ones(1,2),bounds(:,2),[0,0],"LineWidth",1,"Color","k")

% faces with constant x
p(5) = patch(bounds(1,1)*ones(4,1),[bounds(:,2);flip(bounds(:,2))],repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
p(6) = patch(bounds(2,1)*ones(4,1),[bounds(:,2);flip(bounds(:,2))],repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
plot3(bounds(:,1),[0,0],[0,0],"LineWidth",1,"Color","k")

xlim([0 3])
ylim([0 10])
zlim([0 2000])
% legend(ax,ll,"FontSize",20)
legend(ax,"off")
grid on
ax.FontSize = 8;

saveFigures(f,save_fig_opts)

