clearvars;

cohort_name = "cohort_221021113614089";
load(sprintf("../data/%s/%s.mat",cohort_name,cohort_name),"ids","lattice_parameters")


f = gobjects(0,1);
nsamps = size(ids,ndims(ids));
tracked = loadTracked(ids);

%% colors and plotting stuff
therapy_colors = lines(4);
% tumor_type_colors = reshape(jet(4),2,2,3);
% tumor_type_colors = reshape([177,90,95;56,49,50;240,207,75;219,207,186]/255,2,2,3);
% tumor_type_colors = reshape([26,35,73;204,104,35;139,61,30;56,99,105]/255,2,2,3);
tumor_type_colors = reshape([131,111,220;223,89,100;246,218,57;186,224,55]/255,2,2,3);

therapy_titles = ["Together","aPD1 First","aFGFR3 First","No Therapy"];
therapy_legend_title = "Therapy Scheduling";
fgfr3_assumption_titles = ["No Effect","Cytotoxic Effect";"Recruitment Effect","Both Effects"];
tumor_type_labels = ["LA","HA";"LA Mutants","HA Mutants"];
tumor_type_legend_title = "Tumor Type";

%%
t = arrayifyNonuniform(tracked,"t");
max_t = max(t,[],'all','omitnan');
if abs(round(max_t)-max_t) < 1e-2
    max_t = round(max_t);
end
NT = arrayifyNonuniform(tracked,"NT");
NI = arrayifyNonuniform(tracked,"NI");
TT = arrayifyNonuniform(tracked,"tumor_types");
ImmClearance = arrayifyNonuniform(tracked,"imm_cleared");

%% tumor size

f(end+1,1) = figure("Name","tumor_size");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(NT(ri,ci,ti,:,:))';
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},max(0,pc{2}),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = "Tumor Volume (#)";
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% la mut size

f(end+1,1) = figure("Name","tumor_size_la_mut");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(TT(ri,ci,ti,:,:,2,1))';
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},max(0,pc{2}),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["Low Antigen","Mutants (#)"];
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% la mut proportion

f(end+1,1) = figure("Name","tumor_prop_la_mut");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(TT(ri,ci,ti,:,:,2,1))';
            nt_temp = nt_temp./(nt_temp + squeeze(NT(ri,ci,ti,:,:))');
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},min(1,max(0,pc{2})),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["Low Antigen","Mutants (proportion)"];
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% la nonmut size

f(end+1,1) = figure("Name","tumor_size_la_nonmut");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(TT(ri,ci,ti,:,:,1,1))';
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},max(0,pc{2}),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["Low Antigen","Nonmutants (#)"];
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% la nonmut proportion

f(end+1,1) = figure("Name","tumor_prop_la_nonmut");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(TT(ri,ci,ti,:,:,1,1))';
            nt_temp = nt_temp./(nt_temp + squeeze(NT(ri,ci,ti,:,:))');
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},min(1,max(0,pc{2})),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["Low Antigen","Nonmutants (proportion)"];
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% ha mut size

f(end+1,1) = figure("Name","tumor_size_ha_mut");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(TT(ri,ci,ti,:,:,2,2))';
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},max(0,pc{2}),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["High Antigen","Mutants (#)"];
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% ha mut proportion

f(end+1,1) = figure("Name","tumor_prop_ha_mut");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(TT(ri,ci,ti,:,:,2,2))';
            nt_temp = nt_temp./(nt_temp + squeeze(NT(ri,ci,ti,:,:))');
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},min(1,max(0,pc{2})),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["High Antigen","Mutants (proportion)"];
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% ha nonmut size

f(end+1,1) = figure("Name","tumor_size_ha_nonmut");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(TT(ri,ci,ti,:,:,1,2))';
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},max(0,pc{2}),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["High Antigen","Nonmutants (#)"];
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% ha nonmut proportion

f(end+1,1) = figure("Name","tumor_prop_ha_nonmut");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(TT(ri,ci,ti,:,:,1,2))';
            nt_temp = nt_temp./(nt_temp + squeeze(NT(ri,ci,ti,:,:))');
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},min(1,max(0,pc{2})),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["High Antigen","Nonmutants (proportion)"];
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% ctl size

f(end+1,1) = figure("Name","ctl_size");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(NI(ri,ci,ti,:,:))';
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},max(0,pc{2}),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = "CTL Infiltrate (#)";
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% ctl proportion

f(end+1,1) = figure("Name","ctl_proportion");
ax = gobjects(2,2);
pp = gobjects(2,2,4);
ll = gobjects(2,2,4);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";

        for ti = 1:4 % therapy index
            t_temp = squeeze(t(ri,ci,ti,:,:))';
            nt_temp = squeeze(NI(ri,ci,ti,:,:))';
            nt_temp = nt_temp./(nt_temp + squeeze(NT(ri,ci,ti,:,:))');
            [xx,yy,pc] = my_patchPlot(t_temp,nt_temp,true);
            pp(ri,ci,ti) = patch(ax(ri,ci),pc{1},min(1,max(0,pc{2})),therapy_colors(ti,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",therapy_titles(ti));
            ll(ri,ci,ti) = plot(ax(ri,ci),xx,yy,"Color",therapy_colors(ti,:),"LineWidth",1.25,"DisplayName",therapy_titles(ti));

        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["CTL Infiltrate","(proportion of all cells)"];
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t],"YLim",[0 0.2])
normalizeYLims(ax)
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = therapy_legend_title;
f(end).Position(3:4) = [735 420];

%% tumor type immune clearance, no therapy

f(end+1,1) = figure("Name","tumor_imm_clearance_no_therapy");
ax = gobjects(2,2);
pp = gobjects(2,2,2,2);
ll = gobjects(2,2,2,2);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";
        t_temp = squeeze(t(ri,ci,1,:,:))';
        for mut_ind = 1:2
            for ant_ind = 1:2
                nt_temp = squeeze(ImmClearance(ri,ci,4,:,2:end,mut_ind,ant_ind))';
                nt_temp = nt_temp./(squeeze(TT(ri,ci,4,:,1:end-1,mut_ind,ant_ind))' .* diff(t_temp,1,1));
                for si = nsamps:-1:1
                    r(si)=ksr(squeeze(t_temp(1:end-1,si)),nt_temp(:,si),4/24,ceil(max_t*4));
                end
%                 [xx,yy,pc] = my_patchPlot(t_temp(1:end-1,:),nt_temp,true);
                [xx,yy,pc] = my_patchPlot(arrayify(r,"x",1),arrayify(r,"f",1),true);
                pp(ri,ci,mut_ind,ant_ind) = patch(ax(ri,ci),pc{1},max(0,pc{2}),tumor_type_colors(mut_ind,ant_ind,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",tumor_type_labels(mut_ind,ant_ind));
                ll(ri,ci,mut_ind,ant_ind) = plot(ax(ri,ci),xx,yy,"Color",tumor_type_colors(mut_ind,ant_ind,:),"LineWidth",1.25,"DisplayName",tumor_type_labels(mut_ind,ant_ind));
            end
        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["Clearance by","CTL Rate (d^{-1})"];
        ax(ri,ci).YLim(2) = min(2,ax(ri,ci).YLim(2));
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = tumor_type_legend_title;
f(end).Position(3:4) = [735 420];

%% tumor type immune clearance, pd1 first

f(end+1,1) = figure("Name","tumor_imm_clearance_pd1_first");
ax = gobjects(2,2);
pp = gobjects(2,2,2,2);
ll = gobjects(2,2,2,2);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";
        t_temp = squeeze(t(ri,ci,1,:,:))';
        for mut_ind = 1:2
            for ant_ind = 1:2
                nt_temp = squeeze(ImmClearance(ri,ci,3,:,2:end,mut_ind,ant_ind))';
                nt_temp = nt_temp./(squeeze(TT(ri,ci,3,:,1:end-1,mut_ind,ant_ind))' .* diff(t_temp,1,1));
                for si = nsamps:-1:1
                    r(si)=ksr(squeeze(t_temp(1:end-1,si)),nt_temp(:,si),4/24,ceil(max_t*4));
                end
%                 [xx,yy,pc] = my_patchPlot(t_temp(1:end-1,:),nt_temp,true);
                [xx,yy,pc] = my_patchPlot(arrayify(r,"x",1),arrayify(r,"f",1),true);
                pp(ri,ci,mut_ind,ant_ind) = patch(ax(ri,ci),pc{1},max(0,pc{2}),tumor_type_colors(mut_ind,ant_ind,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",tumor_type_labels(mut_ind,ant_ind));
                ll(ri,ci,mut_ind,ant_ind) = plot(ax(ri,ci),xx,yy,"Color",tumor_type_colors(mut_ind,ant_ind,:),"LineWidth",1.25,"DisplayName",tumor_type_labels(mut_ind,ant_ind));
            end
        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["Clearance by","CTL Rate (d^{-1})"];
        ax(ri,ci).YLim(2) = min(2,ax(ri,ci).YLim(2));
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = tumor_type_legend_title;
f(end).Position(3:4) = [735 420];

%% tumor type immune clearance, fgfr3 first

f(end+1,1) = figure("Name","tumor_imm_clearance_fgfr3_first");
ax = gobjects(2,2);
pp = gobjects(2,2,2,2);
ll = gobjects(2,2,2,2);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";
        t_temp = squeeze(t(ri,ci,1,:,:))';
        for mut_ind = 1:2
            for ant_ind = 1:2
                nt_temp = squeeze(ImmClearance(ri,ci,3,:,2:end,mut_ind,ant_ind))';
                nt_temp = nt_temp./(squeeze(TT(ri,ci,3,:,1:end-1,mut_ind,ant_ind))' .* diff(t_temp,1,1));
                for si = nsamps:-1:1
                    r(si)=ksr(squeeze(t_temp(1:end-1,si)),nt_temp(:,si),4/24,ceil(max_t*4));
                end
%                 [xx,yy,pc] = my_patchPlot(t_temp(1:end-1,:),nt_temp,true);
                [xx,yy,pc] = my_patchPlot(arrayify(r,"x",1),arrayify(r,"f",1),true);
                pp(ri,ci,mut_ind,ant_ind) = patch(ax(ri,ci),pc{1},max(0,pc{2}),tumor_type_colors(mut_ind,ant_ind,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",tumor_type_labels(mut_ind,ant_ind));
                ll(ri,ci,mut_ind,ant_ind) = plot(ax(ri,ci),xx,yy,"Color",tumor_type_colors(mut_ind,ant_ind,:),"LineWidth",1.25,"DisplayName",tumor_type_labels(mut_ind,ant_ind));
            end
        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["Clearance by","CTL Rate (d^{-1})"];
        ax(ri,ci).YLim(2) = min(2,ax(ri,ci).YLim(2));
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = tumor_type_legend_title;
f(end).Position(3:4) = [735 420];

%% tumor type immune clearance, therapy together

f(end+1,1) = figure("Name","tumor_imm_clearance_together");
ax = gobjects(2,2);
pp = gobjects(2,2,2,2);
ll = gobjects(2,2,2,2);

for ri = 1:2
    for ci = 1:2
        ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci]));
        ax(ri,ci).NextPlot = "add";
        t_temp = squeeze(t(ri,ci,1,:,:))';
        for mut_ind = 1:2
            for ant_ind = 1:2
                nt_temp = squeeze(ImmClearance(ri,ci,1,:,2:end,mut_ind,ant_ind))';
                nt_temp = nt_temp./(squeeze(TT(ri,ci,1,:,1:end-1,mut_ind,ant_ind))' .* diff(t_temp,1,1));
                for si = nsamps:-1:1
                    r(si)=ksr(squeeze(t_temp(1:end-1,si)),nt_temp(:,si),4/24,ceil(max_t*4));
                end
%                 [xx,yy,pc] = my_patchPlot(t_temp(1:end-1,:),nt_temp,true);
                [xx,yy,pc] = my_patchPlot(arrayify(r,"x",1),arrayify(r,"f",1),true);
                pp(ri,ci,mut_ind,ant_ind) = patch(ax(ri,ci),pc{1},max(0,pc{2}),tumor_type_colors(mut_ind,ant_ind,:),"FaceAlpha",0.15,"EdgeColor","none","DisplayName",tumor_type_labels(mut_ind,ant_ind));
                ll(ri,ci,mut_ind,ant_ind) = plot(ax(ri,ci),xx,yy,"Color",tumor_type_colors(mut_ind,ant_ind,:),"LineWidth",1.25,"DisplayName",tumor_type_labels(mut_ind,ant_ind));
            end
        end

        title(ax(ri,ci),fgfr3_assumption_titles(ri,ci))
        ax(ri,ci).XLabel.String = "Time (d)";
        ax(ri,ci).YLabel.String = ["Clearance by","CTL Rate (d^{-1})"];
        ax(ri,ci).YLim(2) = min(2,ax(ri,ci).YLim(2));
    end
end
set(ax,"FontSize",16,"XLim",[0 max_t])
L = legend(ax(1,1),reshape(ll(1,1,:),[],1),"Location","best");
L.Title.String = tumor_type_legend_title;
f(end).Position(3:4) = [735 420];


%% print
printFigures(f,cohort_name,[],true);


