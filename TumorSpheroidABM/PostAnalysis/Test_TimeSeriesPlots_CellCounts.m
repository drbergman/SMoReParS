clearvars;
addpath("~/Documents/MATLAB/myfunctions/")

cohort_name = "cohort_221012070831969";
load(sprintf("../data/%s/%s.mat",cohort_name,cohort_name))

%%
f = gobjects(0,1);
fgfr3_effects = ["No effect","R (recruit)","C (cytotoxic)","R+C"];
n_cohorts = numel(fgfr3_effects);
all_nr = ceil(sqrt(n_cohorts));
all_nc = ceil(n_cohorts/all_nr);
tumor_type = ["LA","HA";"LA + Mut","HA + Mut"];


%% all tumor sizes
f(end+1,1) = figure("Name","tumor_size_all");
ax = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        ax(recruit_ind,cytotox_ind) = subplot(2,2,r2c(2,2,[recruit_ind,cytotox_ind])); hold on;
        ax(recruit_ind,cytotox_ind).XLim(2) = eps;
        for si = 1:nsamps_per_condition
            plot(tracked(recruit_ind,cytotox_ind,si).t,tracked(recruit_ind,cytotox_ind,si).NT)
            ax(recruit_ind,cytotox_ind).XLim(2) = max(ax(recruit_ind,cytotox_ind).XLim(2),tracked(recruit_ind,cytotox_ind,si).t(end));
        end
        xlabel(ax(recruit_ind,cytotox_ind),"Time (days)")
        ylabel(ax(recruit_ind,cytotox_ind),["Tumor Size","(cells)"])
        title(ax(recruit_ind,cytotox_ind),fgfr3_effects(recruit_ind + 2*(cytotox_ind-1)))
    end
end

set(ax,"FontSize",16)
normalizeXLims(f(end))
normalizeYLims(f(end))
f(end).Position = [840   861   600   297];

%% patch tumor size
t = arrayifyNonuniform(tracked,"t");
t = permute(t,[4,1,2,3]);
NT = arrayifyNonuniform(tracked,"NT");
NT = permute(NT,[4,1,2,3]);
colors = lines(4);
f(end+1,1) = figure("Name","tumor_size"); hold on;
pp = gobjects(2,2);
ll = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        [xx,yy,pc] = my_patchPlot(squeeze(t(:,recruit_ind,cytotox_ind,:)),squeeze(NT(:,recruit_ind,cytotox_ind,:)),true);
        color_ind = sub2ind([2,2],recruit_ind,cytotox_ind);
        pp(recruit_ind,cytotox_ind) = patch(pc{1},max(0,pc{2}),colors(color_ind,:),'FaceAlpha',0.15,'EdgeColor','black');
        ll(recruit_ind,cytotox_ind) = plot(xx,yy,'Color',colors(color_ind,:),"LineWidth",1.25);
    end
end

set(gca,'FontSize',16)
L = legend(ll(:),fgfr3_effects,"Location","NorthWest");
L.Title.String = ["FGFR3 Effect","on CTLs"];
ylabel("Tumor Size (cells)")
xlabel("Time (days)")
title("FGFR3 Effects on CTLs Change Tumor Growth")
f(end).Position = [768   914   632   244];

%% patch tumor types size
t = arrayifyNonuniform(tracked,"t");
t = reshape(t,n_cohorts,nsamps_per_condition,[]);
TumTypes = arrayifyNonuniform(tracked,"tumor_types");
TumTypes = reshape(TumTypes,n_cohorts,nsamps_per_condition,[],2,2);
colors = lines(n_cohorts);
f(end+1,1) = figure("Name","tumor_size_by_type");
ax = gobjects(2,2);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);
for mut_ind = 1:2
    for ant_ind = 1:2
        ax(mut_ind,ant_ind) = subplot(2,2,r2c(2,2,[mut_ind,ant_ind]));
        hold on;
        for ci = 1:n_cohorts
            [xx,yy,pc] = my_patchPlot(squeeze(t(ci,:,:))',squeeze(TumTypes(ci,:,:,mut_ind,ant_ind))',true);
            pp(ci,mut_ind,ant_ind) = patch(ax(mut_ind,ant_ind),pc{1},max(0,pc{2}),colors(ci,:),'FaceAlpha',0.15,'EdgeColor','black');
            ll(ci,mut_ind,ant_ind) = plot(ax(mut_ind,ant_ind),xx,yy,'Color',colors(ci,:),"LineWidth",1.25);
        end
        title(ax(mut_ind,ant_ind),tumor_type(mut_ind,ant_ind))
    end
end

% normalizeYLims(f(end))
set(ax,'FontSize',16)
L = legend(ax(1,1),ll(:,1,1),fgfr3_effects(:),"Location","best");
L.Title.String = ["FGFR3 Effect","on CTLs"];
ylabel(ax,"Size (#)")
xlabel(ax,"Time (days)")
f(end).Position(3:4) = [1263 516];

%% patch tumor types size by condition
t = arrayifyNonuniform(tracked,"t");
t = reshape(t,n_cohorts,nsamps_per_condition,[]);
max_t = max(t,[],'all');
TumTypes = arrayifyNonuniform(tracked,"tumor_types");
TumTypes = reshape(TumTypes,n_cohorts,nsamps_per_condition,[],2,2);
type_colors = jet(4);
f(end+1,1) = figure("Name","tumor_size_by_type_by_condition");
ax = gobjects(n_cohorts,1);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);
for ci = 1:n_cohorts
    ax(ci) = subplot(all_nr,all_nc,r2c(all_nr,all_nc,ci));
    hold on;
    for mut_ind = 1:2
        for ant_ind = 1:2
            [xx,yy,pc] = my_patchPlot(squeeze(t(ci,:,:))',squeeze(TumTypes(ci,:,:,mut_ind,ant_ind))',true);
            pp(ci,mut_ind,ant_ind) = patch(ax(ci),pc{1},max(0,pc{2}),type_colors(mut_ind+(ant_ind-1)*2,:),'FaceAlpha',0.15,'EdgeColor','black');
            ll(ci,mut_ind,ant_ind) = plot(ax(ci),xx,yy,'Color',type_colors(mut_ind+(ant_ind-1)*2,:),"LineWidth",1.25);
        end
        title(ax(ci),fgfr3_effects(ci))
    end
end

% normalizeYLims(f(end))
set(ax,"XLim",[0 max_t])
set(ax,'FontSize',16)
L = legend(ax(1,1),reshape(ll(1,:,:),[],1),tumor_type(:),"Location","best");
L.Title.String = ["FGFR3 Effect","on CTLs"];
ylabel(ax,"Size (#)")
xlabel(ax,"Time (days)")
f(end).Position(3:4) = [1263 516];

%% patch tumor types proportion
t = arrayifyNonuniform(tracked,"t");
t = reshape(t,n_cohorts,nsamps_per_condition,[]);
TumTypes = arrayifyNonuniform(tracked,"tumor_types");
TumTypes = reshape(TumTypes,n_cohorts,nsamps_per_condition,[],2,2);
colors = lines(n_cohorts);
f(end+1,1) = figure("Name","tumor_prop_by_type");
ax = gobjects(2,2);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);
for mut_ind = 1:2
    for ant_ind = 1:2
        ax(mut_ind,ant_ind) = subplot(2,2,r2c(2,2,[mut_ind,ant_ind]));
        hold on;
        for ci = 1:n_cohorts
            [xx,yy,pc] = my_patchPlot(squeeze(t(ci,:,:))',squeeze(TumTypes(ci,:,:,mut_ind,ant_ind))'./squeeze(sum(TumTypes(ci,:,:,:,:),4:5))',true);
            pp(ci,mut_ind,ant_ind) = patch(ax(mut_ind,ant_ind),pc{1},min(1,max(0,pc{2})),colors(ci,:),'FaceAlpha',0.15,'EdgeColor','black');
            ll(ci,mut_ind,ant_ind) = plot(ax(mut_ind,ant_ind),xx,yy,'Color',colors(ci,:),"LineWidth",1.25);
        end
        title(ax(mut_ind,ant_ind),tumor_type(mut_ind,ant_ind))
    end
end

normalizeYLims(f(end))
set(ax,'FontSize',16)
L = legend(ax(1,1),ll(:,1,1),fgfr3_effects(:),"Location","best");
L.Title.String = ["FGFR3 Effect","on CTLs"];
ylabel(ax,"Proportion")
xlabel(ax,"Time (days)")
f(end).Position(3:4) = [1263 516];

%% patch tumor types proportion by condition
t = arrayifyNonuniform(tracked,"t");
t = reshape(t,n_cohorts,nsamps_per_condition,[]);
max_t = max(t,[],'all');
TumTypes = arrayifyNonuniform(tracked,"tumor_types");
TumTypes = reshape(TumTypes,n_cohorts,nsamps_per_condition,[],2,2);
type_colors = jet(4);
f(end+1,1) = figure("Name","tumor_prop_by_type_by_condition");
ax = gobjects(n_cohorts,1);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);
for ci = 1:n_cohorts
    ax(ci) = subplot(all_nr,all_nc,r2c(all_nr,all_nc,ci));
    hold on;
    for mut_ind = 1:2
        for ant_ind = 1:2
            [xx,yy,pc] = my_patchPlot(squeeze(t(ci,:,:))',squeeze(TumTypes(ci,:,:,mut_ind,ant_ind))'./squeeze(sum(TumTypes(ci,:,:,:,:),4:5))',true);
            pp(ci,mut_ind,ant_ind) = patch(ax(ci),pc{1},min(1,max(0,pc{2})),type_colors(mut_ind+(ant_ind-1)*2,:),'FaceAlpha',0.15,'EdgeColor','black');
            ll(ci,mut_ind,ant_ind) = plot(ax(ci),xx,yy,'Color',type_colors(mut_ind+(ant_ind-1)*2,:),"LineWidth",1.25);
        end
        title(ax(ci),fgfr3_effects(ci))
    end
end

normalizeYLims(f(end))
set(ax,"XLim",[0 max_t])
set(ax,'FontSize',16)
L = legend(ax(1,1),reshape(ll(1,:,:),[],1),tumor_type(:),"Location","best");
L.Title.String = ["FGFR3 Effect","on CTLs"];
ylabel(ax,"Proportion")
xlabel(ax,"Time (days)")
f(end).Position(3:4) = [1263 516];

%% all HA sizes
f(end+1,1) = figure("Name","tumor_ha_size_all");
ax = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        ax(recruit_ind,cytotox_ind) = subplot(2,2,r2c(2,2,[recruit_ind,cytotox_ind])); hold on;
        ax(recruit_ind,cytotox_ind).XLim(2) = eps;
        for si = 1:nsamps_per_condition
            plot(tracked(recruit_ind,cytotox_ind,si).t,sum(tracked(recruit_ind,cytotox_ind,si).tumor_types(:,:,2),2))
            ax(recruit_ind,cytotox_ind).XLim(2) = max(ax(recruit_ind,cytotox_ind).XLim(2),tracked(recruit_ind,cytotox_ind,si).t(end));
        end
        xlabel(ax(recruit_ind,cytotox_ind),"Time (days)")
        ylabel(ax(recruit_ind,cytotox_ind),["HA Size","(cells)"])
        title(ax(recruit_ind,cytotox_ind),fgfr3_effects(recruit_ind + 2*(cytotox_ind-1)))
    end
end

set(ax,"FontSize",16)
normalizeXLims(f(end))
normalizeYLims(f(end))
f(end).Position = [840   861   600   297];

%% all HA proportions
f(end+1,1) = figure("Name","tumor_ha_prop_all");
ax = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        ax(recruit_ind,cytotox_ind) = subplot(2,2,r2c(2,2,[recruit_ind,cytotox_ind])); hold on;
        ax(recruit_ind,cytotox_ind).XLim(2) = eps;
        for si = 1:nsamps_per_condition
            plot(tracked(recruit_ind,cytotox_ind,si).t,sum(tracked(recruit_ind,cytotox_ind,si).tumor_types(:,:,2),2)./tracked(recruit_ind,cytotox_ind,si).NT)
            ax(recruit_ind,cytotox_ind).XLim(2) = max(ax(recruit_ind,cytotox_ind).XLim(2),tracked(recruit_ind,cytotox_ind,si).t(end));
        end
        xlabel(ax(recruit_ind,cytotox_ind),"Time (days)")
        ylabel(ax(recruit_ind,cytotox_ind),"HA Proportion")
        title(ax(recruit_ind,cytotox_ind),fgfr3_effects(recruit_ind + 2*(cytotox_ind-1)))
    end
end

set(ax,"FontSize",16)
normalizeXLims(f(end))
normalizeYLims(f(end))
f(end).Position = [840   861   600   297];

%% all LA sizes
f(end+1,1) = figure("Name","tumor_la_size_all");
ax = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        ax(recruit_ind,cytotox_ind) = subplot(2,2,r2c(2,2,[recruit_ind,cytotox_ind])); hold on;
        ax(recruit_ind,cytotox_ind).XLim(2) = eps;
        for si = 1:nsamps_per_condition
            plot(tracked(recruit_ind,cytotox_ind,si).t,sum(tracked(recruit_ind,cytotox_ind,si).tumor_types(:,:,1),2))
            ax(recruit_ind,cytotox_ind).XLim(2) = max(ax(recruit_ind,cytotox_ind).XLim(2),tracked(recruit_ind,cytotox_ind,si).t(end));
        end
        xlabel(ax(recruit_ind,cytotox_ind),"Time (days)")
        ylabel(ax(recruit_ind,cytotox_ind),["LA Size","(cells)"])
        title(ax(recruit_ind,cytotox_ind),fgfr3_effects(recruit_ind + 2*(cytotox_ind-1)))
    end
end

set(ax,"FontSize",16)
normalizeXLims(f(end))
normalizeYLims(f(end))
f(end).Position = [840   861   600   297];

%% patch HA size
t = arrayifyNonuniform(tracked,"t");
t = permute(t,[4,1,2,3]);
TumorTypes = arrayifyNonuniform(tracked,"tumor_types");
HA = sum(TumorTypes(:,:,:,:,:,2),5);
HA = permute(HA,[4,1,2,3]);
colors = lines(4);
f(end+1,1) = figure("Name","tumor_ha_size"); hold on;
pp = gobjects(2,2);
ll = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        [xx,yy,pc] = my_patchPlot(squeeze(t(:,recruit_ind,cytotox_ind,:)),squeeze(HA(:,recruit_ind,cytotox_ind,:)),true);
        color_ind = sub2ind([2,2],recruit_ind,cytotox_ind);
        pp(recruit_ind,cytotox_ind) = patch(pc{1},max(0,pc{2}),colors(color_ind,:),'FaceAlpha',0.15,'EdgeColor','black');
        ll(recruit_ind,cytotox_ind) = plot(xx,yy,'Color',colors(color_ind,:),"LineWidth",1.25);
    end
end

set(gca,'FontSize',16)
L = legend(ll(:),fgfr3_effects,"Location","NorthWest");
L.Title.String = ["FGFR3 Effect","on CTLs"];
ylabel("HA Size (cells)")
xlabel("Time (days)")
title("FGFR3 Effects on CTLs Change HA Tumor Growth")
f(end).Position = [768   914   650   244];

%% patch HA proportion
t = arrayifyNonuniform(tracked,"t");
t = permute(t,[4,1,2,3]);
TumorTypes = arrayifyNonuniform(tracked,"tumor_types");
HA = sum(TumorTypes(:,:,:,:,:,2),5);
HA = permute(HA,[4,1,2,3]);
NT = arrayifyNonuniform(tracked,"NT");
NT = permute(NT,[4,1,2,3]);
colors = lines(4);
f(end+1,1) = figure("Name","tumor_ha_prop"); hold on;
pp = gobjects(2,2);
ll = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        [xx,yy,pc] = my_patchPlot(squeeze(t(:,recruit_ind,cytotox_ind,:)),squeeze(HA(:,recruit_ind,cytotox_ind,:))./squeeze(NT(:,recruit_ind,cytotox_ind,:)),true);
        color_ind = sub2ind([2,2],recruit_ind,cytotox_ind);
        pp(recruit_ind,cytotox_ind) = patch(pc{1},max(0,pc{2}),colors(color_ind,:),'FaceAlpha',0.15,'EdgeColor','black');
        ll(recruit_ind,cytotox_ind) = plot(xx,yy,'Color',colors(color_ind,:),"LineWidth",1.25);
    end
end

set(gca,'FontSize',16)
L = legend(ll(:),fgfr3_effects,"Location","NorthWest");
L.Title.String = ["FGFR3 Effect","on CTLs"];
L.Position = [0.7162    0.4985    0.1877    0.4067];
ylabel("HA Proportion")
xlabel("Time (days)")
title("FGFR3 Effects on CTLs Change HA Proportion")
f(end).Position = [768   874   650   284];

%% patch LA size
t = arrayifyNonuniform(tracked,"t");
t = permute(t,[4,1,2,3]);
NT = arrayifyNonuniform(tracked,"NT");
TumorTypes = arrayifyNonuniform(tracked,"tumor_types");
HA = sum(TumorTypes(:,:,:,:,:,2),5);
LA = permute(NT-HA,[4,1,2,3]);
colors = lines(4);
f(end+1,1) = figure("Name","tumor_la_size"); hold on;
pp = gobjects(2,2);
ll = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        [xx,yy,pc] = my_patchPlot(squeeze(t(:,recruit_ind,cytotox_ind,:)),squeeze(LA(:,recruit_ind,cytotox_ind,:)),true);
        color_ind = sub2ind([2,2],recruit_ind,cytotox_ind);
        pp(recruit_ind,cytotox_ind) = patch(pc{1},max(0,pc{2}),colors(color_ind,:),'FaceAlpha',0.15,'EdgeColor','black');
        ll(recruit_ind,cytotox_ind) = plot(xx,yy,'Color',colors(color_ind,:),"LineWidth",1.25);
    end
end

set(gca,'FontSize',16)
L = legend(ll(:),fgfr3_effects,"Location","NorthWest");
L.Title.String = ["FGFR3 Effect","on CTLs"];
ylabel("LA Size (cells)")
xlabel("Time (days)")
title("FGFR3 Effects on CTLs Change LA Tumor Growth")
f(end).Position = [768   914   650   244];

%% all ctl sizes
f(end+1,1) = figure;
ax = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        ax(recruit_ind,cytotox_ind) = subplot(2,2,r2c(2,2,[recruit_ind,cytotox_ind])); hold on;
        for si = 1:nsamps_per_condition
            plot(tracked(recruit_ind,cytotox_ind,si).t,tracked(recruit_ind,cytotox_ind,si).NI)
        end
    end
end

normalizeXLims(f(end))
normalizeYLims(f(end))

%% patch ctl size
t = arrayifyNonuniform(tracked,"t");
t = permute(t,[4,1,2,3]);
NI = arrayifyNonuniform(tracked,"NI");
NI = permute(NI,[4,1,2,3]);
colors = lines(4);
f(end+1,1) = figure("Name","ctl_size"); hold on;
pp = gobjects(2,2);
ll = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        [xx,yy,pc] = my_patchPlot(squeeze(t(:,recruit_ind,cytotox_ind,:)),squeeze(NI(:,recruit_ind,cytotox_ind,:)),true);
        color_ind = sub2ind([2,2],recruit_ind,cytotox_ind);
        pp(recruit_ind,cytotox_ind) = patch(pc{1},max(0,pc{2}),colors(color_ind,:),'FaceAlpha',0.15,'EdgeColor','black');
        ll(recruit_ind,cytotox_ind) = plot(xx,yy,'Color',colors(color_ind,:),"LineWidth",1.25);
    end
end

set(gca,'FontSize',16)
L = legend(ll(:),fgfr3_effects,"Location","NorthWest");
L.Title.String = ["FGFR3 Effect","on CTLs"];
ylabel("CTL Inifiltrate (cells)")
xlabel("Time (days)")
title("FGFR3 Effects on CTLs Change CTL Dynamics")
f(end).Position = [768   914   632   244];

%% all ctl proportions

propfn = @(x,y) y./(x+y);
f(end+1,1) = figure;
ax = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        ax(recruit_ind,cytotox_ind) = subplot(2,2,r2c(2,2,[recruit_ind,cytotox_ind])); hold on;
        for si = 1:nsamps_per_condition
            plot(tracked(recruit_ind,cytotox_ind,si).t,propfn(tracked(recruit_ind,cytotox_ind,si).NT,tracked(recruit_ind,cytotox_ind,si).NI))
        end
    end
end

normalizeXLims(f(end))
normalizeYLims(f(end))

%% patch ctl proportions
t = arrayifyNonuniform(tracked,"t");
t = permute(t,[4,1,2,3]);
NT = arrayifyNonuniform(tracked,"NT");
NT = permute(NT,[4,1,2,3]);
NI = arrayifyNonuniform(tracked,"NI");
NI = permute(NI,[4,1,2,3]);
Y = NI./(NI+NT);
colors = lines(4);
f(end+1,1) = figure("Name","ctl_proportions"); hold on;
pp = gobjects(2,2);
ll = gobjects(2,2);
for recruit_ind = 1:2
    for cytotox_ind = 1:2
        [xx,yy,pc] = my_patchPlot(squeeze(t(:,recruit_ind,cytotox_ind,:)),squeeze(Y(:,recruit_ind,cytotox_ind,:)),true);
        color_ind = sub2ind([2,2],recruit_ind,cytotox_ind);
        pp(recruit_ind,cytotox_ind) = patch(pc{1},min(1,max(0,pc{2})),colors(color_ind,:),'FaceAlpha',0.15,'EdgeColor','none');
        ll(recruit_ind,cytotox_ind) = plot(xx,yy,'Color',colors(color_ind,:),"LineWidth",1.25);
    end
end

set(gca,'FontSize',16)
L = legend(ll(:),fgfr3_effects,"Location","NorthWest");
L.Title.String = ["FGFR3 Effect","on CTLs"];
ylabel("CTL Proportion")
xlabel("Time (days)")
title("FGFR3 Effects on CTL Proportion")
f(end).Position = [768   914   632   244];

%% print
printFigures(f,cohort_name);
