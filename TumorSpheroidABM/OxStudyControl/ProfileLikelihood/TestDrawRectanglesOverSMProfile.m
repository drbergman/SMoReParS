clearvars;

load("data/Profiles_SMFromABM_New_clean.mat")

n_sm_pars = size(profiles,1);
profiles = reshape(profiles,n_sm_pars,[]);
n_abm_pars = size(profiles,2);

face_color = [.1 .1 .1];
f = openfig("figures/fig/ParTriplesAll_SMFitToData_New.fig");
ax = gca;
set(ax.Children,"LineWidth",2,"Marker","none","LineStyle","--")
view(43,62)
legend("AutoUpdate","off")

show_one_box_at_a_time = true;
xlim([0 3])
zlim([0 2000])
% zlim([0 1000])
% xlim([0 5])
% ylim([0 10])
hold on
line_style_in_box = "-";
line_style_out_box = ":";
line_color_in_box = "";
line_width_in_box = 2;
for i = 1:3
    xx{i} = ax.Children(i).XData;
    yy{i} = ax.Children(i).YData;
    zz{i} = ax.Children(i).ZData;
    colors(i,:) = ax.Children(i).Color;
    disp_names{i} = ax.Children(i).DisplayName;
    ax.Children(i).LineWidth = 1;
    ax.Children(i).LineStyle = line_style_out_box;
end
for i = 1:3
    %% projections onto each axis
    if i==1
        xx_temp = xx{i};
        yy_temp = yL(1)*ones(1,length(xx_temp));
        zz_temp = zL(1)*ones(1,length(xx_temp));
    elseif i==2
        yy_temp = yy{i};
        xx_temp = xL(2)*ones(1,length(yy_temp));
        zz_temp = zL(1)*ones(1,length(yy_temp));
    else
        zz_temp = zz{i};
        xx_temp = xL(1)*ones(1,length(zz_temp));
        yy_temp = yL(1)*ones(1,length(zz_temp));
    end
    plot3(ax,xx_temp,yy_temp,zz_temp,"Color",colors(i,:),"LineWidth",3,"LineStyle",":");

end
for abm_par_ind = 1:n_abm_pars

    bounds = zeros(2,n_sm_pars);
    for sm_par_ind = 1:n_sm_pars
        bounds(:,sm_par_ind) = profiles{sm_par_ind,abm_par_ind}(sm_par_ind,[1,end]);
    end
    % faces with constant z
    p(1) = patch([bounds(:,1);flip(bounds(:,1))],repelem(bounds(:,2),2,1),bounds(1,3)*ones(4,1),face_color,"FaceAlpha",0.25);
    p(2) = patch([bounds(:,1);flip(bounds(:,1))],repelem(bounds(:,2),2,1),bounds(2,3)*ones(4,1),face_color,"FaceAlpha",0.25);
    
    % faces with constant y
    p(3) = patch([bounds(:,1);flip(bounds(:,1))],bounds(1,2)*ones(4,1),repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
    p(4) = patch([bounds(:,1);flip(bounds(:,1))],bounds(2,2)*ones(4,1),repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
    
    % faces with constant x
    p(5) = patch(bounds(1,1)*ones(4,1),[bounds(:,2);flip(bounds(:,2))],repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
    p(6) = patch(bounds(2,1)*ones(4,1),[bounds(:,2);flip(bounds(:,2))],repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);

    drawnow
    if show_one_box_at_a_time
        delete(p)
    end
end
