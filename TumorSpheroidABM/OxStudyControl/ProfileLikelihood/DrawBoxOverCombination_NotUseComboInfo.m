clearvars;

load("data/Profiles_SMFromABM_New_clean.mat")

save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig";"png"];
save_fig_opts.resolution = '-r1200';

n_sm_pars = size(profiles,1);
profiles = reshape(profiles,n_sm_pars,[]);
n_abm_pars = size(profiles,2);

line_style_in_box = "-";
line_style_out_box = ":";
line_color_in_box = "";
line_width_in_box = 1;
%% axis bounds
xL = [0 2];
yL = [2 10];
zL = [0 2000];

%% draw box
face_color = [.1 .1 .1];
f = openfig("figures/fig/ParTriplesAll_SMFitToData_New.fig");
f.Name = "DifferntComboInfoBox";
ax = gca;
set(ax.Children,"LineWidth",0.5,"Marker","none","LineStyle","--")
view(27.6379,12.2458)
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
plot3(xL(1)+[0,0],yL(1)+[0,0],bounds(:,3),"LineWidth",2,"Color","k")

% faces with constant y
p(3) = patch([bounds(:,1);flip(bounds(:,1))],bounds(1,2)*ones(4,1),repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
p(4) = patch([bounds(:,1);flip(bounds(:,1))],bounds(2,2)*ones(4,1),repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
plot3(xL(2)+[0,0],bounds(:,2),zL(1)+[0,0],"LineWidth",2,"Color","k")

% faces with constant x
p(5) = patch(bounds(1,1)*ones(4,1),[bounds(:,2);flip(bounds(:,2))],repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
p(6) = patch(bounds(2,1)*ones(4,1),[bounds(:,2);flip(bounds(:,2))],repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
plot3(bounds(:,1),yL(1)+[0,0],zL(1)+[0,0],"LineWidth",2,"Color","k")

patch(ax,[bounds(:,1);flip(bounds(:,1))],repelem(bounds(:,2),2,1),zL(1)*ones(4,1),face_color,"FaceAlpha",0.25);
patch(ax,[bounds(:,1);flip(bounds(:,1))],yL(2)*ones(4,1),repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);
patch(ax,xL(1)*ones(4,1),[bounds(:,2);flip(bounds(:,2))],repelem(bounds(:,3),2,1),face_color,"FaceAlpha",0.25);

for i = 1:2
    plot3(xx{i},yy{i},zz{i},"Color",colors(i,:),"LineWidth",1)
    % P1 = interp1(yy{i},[xx{i}',yy{i}',zz{i}'],bounds(2,2));
    % I1 = yy{i}>=bounds(2,2);
    % c{1} = [xx{i}(I1),P1(1)];
    % c{2} = [yy{i}(I1),P1(2)];
    % c{3} = [zz{i}(I1),P1(3)];
    % [~,order] = sort(c{i});
    % plot3(c{1}(order),c{2}(order),c{3}(order),"LineStyle",line_style_out_box,"LineWidth",0.5,"Color",colors(i,:));
    % 
    % P2 = interp1(xx{i},[xx{i}',yy{i}',zz{i}'],bounds(2,1));
    % I3 = xx{i}>=bounds(2,1);
    % c{1} = [xx{i}(I3),P2(1)];
    % c{2} = [yy{i}(I3),P2(2)];
    % c{3} = [zz{i}(I3),P2(3)];
    % [~,order] = sort(c{i});
    % plot3(c{1}(order),c{2}(order),c{3}(order),"LineStyle",line_style_out_box,"LineWidth",0.5,"Color",colors(i,:));
    % 
    % I2 = ~I1 & ~I3;
    % c{1} = [xx{i}(I2),P1(1),P2(1)];
    % c{2} = [yy{i}(I2),P1(2),P2(2)];
    % c{3} = [zz{i}(I2),P1(3),P2(3)];
    % [~,order] = sort(c{i});
    % if line_color_in_box=="white"
    %     temp = line_color_in_box;
    % else
    %     temp = colors(i,:);
    % end
    % ll(i) = plot3(c{1}(order),c{2}(order),c{3}(order),"LineStyle",line_style_in_box,"LineWidth",line_width_in_box,"Color",temp,"DisplayName",disp_names{i});

    %% projections onto planes
    other_coord = setdiff(1:2,i);
    for j = 1:numel(other_coord)
        i1 = i;
        i2 = other_coord(j);
        if i1>i2
            temp = i1;
            i1 = i2;
            i2 = temp;
        end
        path_ = [xx{i}',yy{i}',zz{i}'];
        path_ = path_(:,[i1,i2]);
        % inside_box = (path_(:,1) <= bounds(2,i1)) & (path_(:,1) >= bounds(1,i1)) & (path_(:,2) <= bounds(2,i2)) & (path_(:,2) >= bounds(1,i2));
        if i1==1
            xx_temp = path_(:,1);
            if i2==2
                yy_temp = path_(:,2);
                zz_temp = zL(1)*ones(size(path_,1),1);
            else
                zz_temp = path_(:,2);
                yy_temp = yL(2)*ones(size(path_,1),1);
            end
        else
            xx_temp = xL(1)*ones(size(path_,1),1);
            yy_temp = path_(:,1);
            zz_temp = path_(:,2);
        end
        plot3(ax,xx_temp,yy_temp,zz_temp,"Color",colors(i,:),"LineWidth",2,"LineStyle","-.");
    end

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
        xx_temp = xL(2)*ones(1,length(zz_temp));
        yy_temp = yL(1)*ones(1,length(zz_temp));
    end
    plot3(ax,xx_temp,yy_temp,zz_temp,"Color",colors(i,:),"LineWidth",2,"LineStyle",":");

end

f.Position(3) = 3;

xlim(xL)
ylim(yL)
zlim(zL)
% legend(ax,ll,"FontSize",20)
legend(ax,"off")
grid on
ax.FontSize = 8;

%%
margin = struct("left",0.15,"right",.07,"top",0,"bottom",0.14);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);
%%
saveFigures(f,save_fig_opts)
