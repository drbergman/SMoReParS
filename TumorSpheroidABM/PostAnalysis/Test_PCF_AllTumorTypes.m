clearvars;
addpath("~/Documents/MATLAB/myfunctions/")

cohort_name = "cohort_349087890609125";
load(sprintf("../data/%s/pcf_low_antigen_nonmut_to_immune.mat",cohort_name))
out_la_nm = out;
load(sprintf("../data/%s/pcf_high_antigen_nonmut_to_immune.mat",cohort_name))
out_ha_nm = out;
load(sprintf("../data/%s/pcf_low_antigen_mut_to_immune.mat",cohort_name))
out_la_m = out;
load(sprintf("../data/%s/pcf_high_antigen_mut_to_immune.mat",cohort_name))
out_ha_m = out;
load(sprintf("../data/%s/%s.mat",cohort_name,cohort_name))
clear out
f = gobjects(0,1);

%%
tracked = reshape(tracked,[],nsamps_per_condition);
out_la_nm = reshape(out_la_nm,[],nsamps_per_condition);
out_ha_nm = reshape(out_ha_nm,[],nsamps_per_condition);
out_la_m = reshape(out_la_m,[],nsamps_per_condition);
out_ha_m = reshape(out_ha_m,[],nsamps_per_condition);
out_all = {out_la_nm,out_ha_nm;out_la_m,out_ha_m};

n_cohorts = size(tracked,1);

%% all time series plots by assumption
nr = ceil(sqrt(nsamps_per_condition));
nc = ceil(nsamps_per_condition/nr);
ax = gobjects(n_cohorts,nsamps_per_condition);
for ci = 1:n_cohorts
    f(end+1,1) = figure;
    for si = 1:nsamps_per_condition
        ax(ci,si) = subplot(nr,nc,si);
        pcf = out_ha_nm(ci,si).avg'./out_la_nm(ci,si).avg';
        pcf(out_la_nm(ci,si).avg'==0) = NaN;
        pcfTimeSeriesPlot(options.rr,out_la_nm(ci,si).t,pcf)
        drawnow
    end
end

%% all time series plots by sample
nr = 2;
nc = 2;
ax = gobjects(n_cohorts,nsamps_per_condition);
for si = 1:nsamps_per_condition
    f(end+1) = figure;
    for ci = 1:n_cohorts
        ax(ci,si) = subplot(nr,nc,r2c(nr,nc,ci));
        pcf = out_ha_nm(ci,si).avg'./out_la_nm(ci,si).avg';
        pcf(out_la_nm(ci,si).avg'==0) = NaN;
        pcfTimeSeriesPlot(options.rr,out_la_nm(ci,si).t,pcf)
        drawnow
    end
end

%% average time series
nr = 2;
nc = 2;
ax = gobjects(2,2);
f(end+1) = figure("Name","pcf_average_ha2la_nonmuts");
threshold = .0*nsamps_per_condition;
pcf_limits = [Inf,-Inf];
fgfr3_effects = ["No effect","Cytotoxic effect";"Recruit effect","Both effects"];
for ci = 1:n_cohorts
    ax(ci) = subplot(nr,nc,r2c(nr,nc,ci));
    total_avg = zeros(length(options.rr),0);
    total_avg_la = zeros(length(options.rr),0);
    weights = zeros(1,0);
    t = [];
    for si = 1:nsamps_per_condition
        t = unique([t,out_la_nm(ci,si).t]);
        szs = [size(total_avg,2),size(out_la_nm(ci,si).avg,2)]; % number of time points for each so far
        if szs(1)<szs(2)
            weights = weights+1;
            total_avg = ((si-1)*total_avg + out_ha_nm(ci,si).avg(:,1:szs(1)))./weights;
            total_avg_la = ((si-1)*total_avg_la + out_la_nm(ci,si).avg(:,1:szs(1)))./weights;
            weights(end+1:szs(2)) = 1;
            total_avg(:,end+1:szs(2)) = out_ha_nm(ci,si).avg(:,szs(1)+1:end);
            total_avg_la(:,end+1:szs(2)) = out_la_nm(ci,si).avg(:,szs(1)+1:end);
        elseif szs(1)==szs(2)
            weights = weights+1;
            total_avg = ((si-1)*total_avg + out_ha_nm(ci,si).avg)./weights;
            total_avg_la = ((si-1)*total_avg_la + out_la_nm(ci,si).avg)./weights;
        else
            weights(1:szs(2)) = weights(1:szs(2))+1;
            total_avg(:,1:szs(2)) = ((si-1)*total_avg(:,1:szs(2)) + out_ha_nm(ci,si).avg)./weights(1:szs(2));
            total_avg_la(:,1:szs(2)) = ((si-1)*total_avg_la(:,1:szs(2)) + out_la_nm(ci,si).avg)./weights(1:szs(2));
        end
    end
    pcf = total_avg(:,weights>threshold)'./total_avg_la(:,weights>threshold)';
    pcf(total_avg_la(:,weights>threshold)'==0) = NaN;
    pcf_limits(1) = min([pcf_limits(1);pcf(:)],[],'omitnan');
    pcf_limits(2) = max([pcf_limits(2);pcf(:)],[],'omitnan');
    pcfTimeSeriesPlot(options.rr,t(weights>threshold),pcf)
    title(fgfr3_effects(ci))
end

normalizeYLims(f(end))

pcf_limits(2) = min(pcf_limits(2),4);
one_color = 1*ones(1,3);
ngrid = 1001;
max_colors = 101; % max colors for each of positive and negative directions
positive_corr_n = round(ngrid*(pcf_limits(2)-1)/diff(pcf_limits));
negative_corr_n = ngrid-positive_corr_n;

% unique_pos_colors = floor((max_colors+1)*linspace(0,1,positive_corr_n))/(max_colors+1);
% unique_neg_colors = floor((max_colors+1)*linspace(0,1,negative_corr_n))/(max_colors+1);
unique_pos_color_green = floor((max_colors+1)*linspace(0,one_color(2),positive_corr_n))/(max_colors+1);
unique_pos_color_blue = floor((max_colors+1)*linspace(0,one_color(3),positive_corr_n))/(max_colors+1);
unique_neg_colors_red = floor((max_colors+1)*linspace(0,one_color(1),negative_corr_n))/(max_colors+1);
unique_neg_colors_green = floor((max_colors+1)*linspace(0,one_color(2),negative_corr_n))/(max_colors+1);

% positive_colors = [ones(positive_corr_n,1),flip(unique_pos_colors')*[1,1]];
% negative_colors = [unique_neg_colors'*[1,1],ones(negative_corr_n,1)];
positive_colors = [linspace(1,one_color(1),positive_corr_n)',flip(unique_pos_color_green'),flip(unique_pos_color_blue')];
negative_colors = [unique_neg_colors_red',unique_neg_colors_green',linspace(1,one_color(3),negative_corr_n)'];

cmap = [negative_colors;positive_colors];
for axi = 1:numel(ax)
    colormap(ax(axi),cmap);
    clim(ax(axi),pcf_limits)
    c = colorbar(ax(axi));
    c.Label.String = "g_{ha}(r) / g_{la}(r)";
    c.Label.Rotation = 270;
    c.Label.VerticalAlignment = "baseline";
end
set(ax,'FontSize',16)
f(end).Position = [1231         738         872         420];

% %% plot single sample from each
% nr=2;
% nc=2;
% ax = gobjects(2,2);
% f(end+1) = figure;
% tvals = [2,6,15,25,30,50];
% colors = parula(length(tvals));
% for recruit_ind = 1:2
%     for cytotox_ind = 1:2
%         ax(ci) = subplot(nr,nc,r2c(nr,nc,[ci]));
%         for ti = 1:length(tvals)
%             i = find(out(ci,1).t<=tvals(ti),1,'last');
%             xx = [options.rr,flip(options.rr)]';
%             yy = [out(ci,1).avg(:,i);flip(out(ci,1).avg(:,i))] + [out(ci,1).std(:,i);-flip(out(ci,1).std(:,i))];
%             
%             bad_ind = isnan(xx) | isnan(yy);
%             xx(bad_ind) = [];
%             yy(bad_ind) = [];
%             
%             if isempty(xx)
%                 continue;
%             end
%             pp(ci,ti) = patch(xx,max(0,yy),colors(ti,:),'FaceAlpha',0.2);
%         end
% 
% 
%     end
% end

%% average time series by type
nr = 2;
nc = 2;
ax = gobjects(n_cohorts,2,2);
fgfr3_effects = ["No effect","Cytotoxic effect";"Recruit effect","Both effects"];
for ci = 1:n_cohorts
    f(end+1) = figure("Name",sprintf("pcf_average_%s",fgfr3_effects(ci)));
    threshold = .0*nsamps_per_condition;
    pcf_limits = [Inf,-Inf];
    for mut_ind = 1:2
        for ant_ind = 1:2
            ax(ci,mut_ind,ant_ind) = subplot(nr,nc,r2c(nr,nc,[mut_ind,ant_ind]));
            total_avg = zeros(length(options.rr),0);
            weights = zeros(1,0);
            t = [];
            for si = 1:nsamps_per_condition
                t = unique([t,out_all{mut_ind,ant_ind}(ci,si).t]);
                szs = [size(total_avg,2),size(out_all{mut_ind,ant_ind}(ci,si).avg,2)]; % number of time points for each so far
                if szs(1)<szs(2)
                    weights = weights+1;
                    total_avg = ((si-1)*total_avg + out_ha_nm(ci,si).avg(:,1:szs(1)))./weights;
                    total_avg = ((si-1)*total_avg + out_all{mut_ind,ant_ind}(ci,si).avg(:,1:szs(1)))./weights;
                    weights(end+1:szs(2)) = 1;
                    total_avg(:,end+1:szs(2)) = out_ha_nm(ci,si).avg(:,szs(1)+1:end);
                    total_avg(:,end+1:szs(2)) = out_all{mut_ind,ant_ind}(ci,si).avg(:,szs(1)+1:end);
                elseif szs(1)==szs(2)
                    weights = weights+1;
                    total_avg = ((si-1)*total_avg + out_ha_nm(ci,si).avg)./weights;
                    total_avg = ((si-1)*total_avg + out_all{mut_ind,ant_ind}(ci,si).avg)./weights;
                else
                    weights(1:szs(2)) = weights(1:szs(2))+1;
                    total_avg(:,1:szs(2)) = ((si-1)*total_avg(:,1:szs(2)) + out_ha_nm(ci,si).avg)./weights(1:szs(2));
                    total_avg(:,1:szs(2)) = ((si-1)*total_avg(:,1:szs(2)) + out_all{mut_ind,ant_ind}(ci,si).avg)./weights(1:szs(2));
                end
            end
            pcf = total_avg(:,weights>threshold)'./total_avg(:,weights>threshold)';
            pcf(total_avg(:,weights>threshold)'==0) = NaN;
            pcf_limits(1) = min([pcf_limits(1);pcf(:)],[],'omitnan');
            pcf_limits(2) = max([pcf_limits(2);pcf(:)],[],'omitnan');
            pcfTimeSeriesPlot(options.rr,t(weights>threshold),pcf)
            title(fgfr3_effects(ci))
        end
    end

    normalizeYLims(f(end))

    pcf_limits(2) = min(pcf_limits(2),4);
    one_color = 1*ones(1,3);
    ngrid = 1001;
    max_colors = 101; % max colors for each of positive and negative directions
    positive_corr_n = round(ngrid*(pcf_limits(2)-1)/diff(pcf_limits));
    negative_corr_n = ngrid-positive_corr_n;

    % unique_pos_colors = floor((max_colors+1)*linspace(0,1,positive_corr_n))/(max_colors+1);
    % unique_neg_colors = floor((max_colors+1)*linspace(0,1,negative_corr_n))/(max_colors+1);
    unique_pos_color_green = floor((max_colors+1)*linspace(0,one_color(2),positive_corr_n))/(max_colors+1);
    unique_pos_color_blue = floor((max_colors+1)*linspace(0,one_color(3),positive_corr_n))/(max_colors+1);
    unique_neg_colors_red = floor((max_colors+1)*linspace(0,one_color(1),negative_corr_n))/(max_colors+1);
    unique_neg_colors_green = floor((max_colors+1)*linspace(0,one_color(2),negative_corr_n))/(max_colors+1);

    % positive_colors = [ones(positive_corr_n,1),flip(unique_pos_colors')*[1,1]];
    % negative_colors = [unique_neg_colors'*[1,1],ones(negative_corr_n,1)];
    positive_colors = [linspace(1,one_color(1),positive_corr_n)',flip(unique_pos_color_green'),flip(unique_pos_color_blue')];
    negative_colors = [unique_neg_colors_red',unique_neg_colors_green',linspace(1,one_color(3),negative_corr_n)'];

    cmap = [negative_colors;positive_colors];
    for axi = 1:numel(ax)
        colormap(ax(axi),cmap);
        clim(ax(axi),pcf_limits)
        c = colorbar(ax(axi));
        c.Label.String = "g_{ha}(r) / g_{la}(r)";
        c.Label.Rotation = 270;
        c.Label.VerticalAlignment = "baseline";
    end
    set(ax,'FontSize',16)
    f(end).Position = [1231         738         872         420];
end
% %% plot single sample from each
% nr=2;
% nc=2;
% ax = gobjects(2,2);
% f(end+1) = figure;
% tvals = [2,6,15,25,30,50];
% colors = parula(length(tvals));
% for recruit_ind = 1:2
%     for cytotox_ind = 1:2
%         ax(ci) = subplot(nr,nc,r2c(nr,nc,[ci]));
%         for ti = 1:length(tvals)
%             i = find(out(ci,1).t<=tvals(ti),1,'last');
%             xx = [options.rr,flip(options.rr)]';
%             yy = [out(ci,1).avg(:,i);flip(out(ci,1).avg(:,i))] + [out(ci,1).std(:,i);-flip(out(ci,1).std(:,i))];
%             
%             bad_ind = isnan(xx) | isnan(yy);
%             xx(bad_ind) = [];
%             yy(bad_ind) = [];
%             
%             if isempty(xx)
%                 continue;
%             end
%             pp(ci,ti) = patch(xx,max(0,yy),colors(ti,:),'FaceAlpha',0.2);
%         end
% 
% 
%     end
% end

%% arg max patch plots
f(end+1) = figure("Name","pcf_argmax_tumtypes_to_immune");
ax = gobjects(2,2);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);
colors = lines(4);
tumor_type = ["LA","HA";"LA + Mut","HA + Mut"];
fgfr3_effects = ["No effect","C (cytotoxic)";"R (recruit)","R+C"];

for mut_ind = 1:2
    for ant_ind = 1:2
        ax(mut_ind,ant_ind) = subplot(2,2,r2c(2,2,[mut_ind,ant_ind]));
        hold on;
        t = arrayifyNonuniform(out_all{mut_ind,ant_ind},"t");
        max_t = max(t,[],'all');
        for i = 1:numel(out_all{mut_ind,ant_ind})
            [~,out_all{mut_ind,ant_ind}(i).argmax] = max(out_all{mut_ind,ant_ind}(i).avg,[],1);
        end

        argmax = arrayifyNonuniform(out_all{mut_ind,ant_ind},"argmax");
        rmax = NaN(size(argmax));
        for i = 1:numel(argmax)
            if ~isnan(argmax(i))
                rmax(i) = options.rr(argmax(i));
            end
        end
        for ci = 1:n_cohorts
            [xx,yy,pc] = my_patchPlot(squeeze(t(ci,:,:))',squeeze(rmax(ci,:,:))',true);
            color_ind = sub2ind([2,2],ci);
            pp(ci,mut_ind,ant_ind) = patch(ax(mut_ind,ant_ind),pc{1},max(0,pc{2}),colors(color_ind,:),'FaceAlpha',0.15,'EdgeColor','none','DisplayName',fgfr3_effects(ci));
            ll(ci,mut_ind,ant_ind) = plot(ax(mut_ind,ant_ind),xx,yy,'Color',colors(color_ind,:),"LineWidth",1.25,'DisplayName',fgfr3_effects(ci));
        end
        title(tumor_type(mut_ind,ant_ind))
        xlim(ax(mut_ind,ant_ind),[0 max_t])
    end
end
normalizeXLims(ax)
xlabel(ax,'Time (days)')
ylabel(ax,'r (cell widths)')
% title(["Radius of Maximal Clustering","of CTLs near HA Cells"])
set(ax,'FontSize',16)
L = legend(reshape(ll(:,1,1),[],1),"Location","best");
L.Title.String = ["FGFR3 Effect","on CTLs"];
L.Position = [0.1570    0.6845    0.0940    0.2310];
f(end).Position(3:4) = [1298 500];

%% arg max patch plots by condition
f(end+1) = figure("Name","pcf_argmax_tumtypes_to_immune_by_condition");
ax = gobjects(n_cohorts,1);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);
type_colors = jet(4);
out_all = {out_la_nm,out_ha_nm;out_la_m,out_ha_m};
tumor_type = ["LA","HA";"LA + Mut","HA + Mut"];
fgfr3_effects = ["No effect","C (cytotoxic)";"R (recruit)","R+C"];

all_nr = ceil(sqrt(n_cohorts));
all_nc = ceil(n_cohorts/all_nr);
for ci = 1:n_cohorts
    ax(ci) = subplot(all_nr,all_nc,r2c(all_nr,all_nc,ci));
    hold on
    title(ax(ci),fgfr3_effects(ci))
end

for mut_ind = 1:2
    for ant_ind = 1:2
        t = arrayifyNonuniform(out_all{mut_ind,ant_ind},"t");
        max_t = max(t,[],'all');
        for i = 1:numel(out_all{mut_ind,ant_ind})
            [~,out_all{mut_ind,ant_ind}(i).argmax] = max(out_all{mut_ind,ant_ind}(i).avg,[],1);
        end
        argmax = arrayifyNonuniform(out_all{mut_ind,ant_ind},"argmax");
        rmax = NaN(size(argmax));
        for i = 1:numel(argmax)
            if ~isnan(argmax(i))
                rmax(i) = options.rr(argmax(i));
            end
        end
        for ci = 1:n_cohorts
            [xx,yy,pc] = my_patchPlot(squeeze(t(ci,:,:))',squeeze(rmax(ci,:,:))',true);
            color_ind = sub2ind([2,2],mut_ind,ant_ind);
            pp(ci,mut_ind,ant_ind) = patch(ax(ci),pc{1},max(0,pc{2}),type_colors(color_ind,:),'FaceAlpha',0.15,'EdgeColor','none','DisplayName',tumor_type(mut_ind,ant_ind));
            ll(ci,mut_ind,ant_ind) = plot(ax(ci),xx,yy,'Color',type_colors(color_ind,:),"LineWidth",1.25,'DisplayName',tumor_type(mut_ind,ant_ind));
            xlim(ax(ci),[0 max_t])
        end
    end
end
normalizeXLims(ax)
xlabel(ax,'Time (days)')
ylabel(ax,'r (cell widths)')
% title(["Radius of Maximal Clustering","of CTLs near HA Cells"])
set(ax,'FontSize',16)
L = legend(reshape(ll(1,:,:),[],1),"Location","best");
L.Title.String = "Tumor Type";

%% r in expectation patch plots
f(end+1) = figure("Name","pcf_rint_tumtypes_to_immune");
ax = gobjects(2,2);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);
colors = lines(4);
out_all = {out_la_nm,out_ha_nm;out_la_m,out_ha_m};
tumor_type = ["LA","HA";"LA + Mut","HA + Mut"];
fgfr3_effects = ["No effect","C (cytotoxic)";"R (recruit)","R+C"];

for mut_ind = 1:2
    for ant_ind = 1:2
        ax(mut_ind,ant_ind) = subplot(2,2,r2c(2,2,[mut_ind,ant_ind]));
        hold on;
        t = arrayifyNonuniform(out_all{mut_ind,ant_ind},"t");
        max_t = max(t,[],'all');
        for i = 1:numel(out_all{mut_ind,ant_ind})
            out_all{mut_ind,ant_ind}(i).rint = sum(out_all{mut_ind,ant_ind}(i).avg.*(options.rr.^1)',1,'omitnan')./sum(out_all{mut_ind,ant_ind}(i).avg,1,'omitnan');
        end

        rint = arrayifyNonuniform(out_all{mut_ind,ant_ind},"rint");
        for ci = 1:n_cohorts
            [xx,yy,pc] = my_patchPlot(squeeze(t(ci,:,:))',squeeze(rint(ci,:,:))',true);
            color_ind = sub2ind([2,2],ci);
            pp(ci,mut_ind,ant_ind) = patch(ax(mut_ind,ant_ind),pc{1},max(0,pc{2}),colors(color_ind,:),'FaceAlpha',0.15,'EdgeColor','none','DisplayName',fgfr3_effects(ci));
            ll(ci,mut_ind,ant_ind) = plot(ax(mut_ind,ant_ind),xx,yy,'Color',colors(color_ind,:),"LineWidth",1.25,'DisplayName',fgfr3_effects(ci));
        end
        title(tumor_type(mut_ind,ant_ind))
        xlim(ax(mut_ind,ant_ind),[0 max_t])
    end
end
normalizeXLims(ax)
xlabel(ax,'Time (days)')
ylabel(ax,'r (cell widths)')
% title(["Radius of Maximal Clustering","of CTLs near HA Cells"])
set(ax,'FontSize',16)
L = legend(reshape(ll(:,1,1),[],1),"Location","best");
L.Title.String = ["FGFR3 Effect","on CTLs"];
L.Position = [0.1632    0.7185    0.0940    0.2091];
f(end).Position(3:4) = [1298 706];

%% r in expectation patch plots by condition
f(end+1) = figure("Name","pcf_rint_tumtypes_to_immune_by_condition");
ax = gobjects(n_cohorts,1);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);
type_colors = jet(4);
out_all = {out_la_nm,out_ha_nm;out_la_m,out_ha_m};
tumor_type = ["LA","HA";"LA + Mut","HA + Mut"];
fgfr3_effects = ["No effect","C (cytotoxic)";"R (recruit)","R+C"];

all_nr = ceil(sqrt(n_cohorts));
all_nc = ceil(n_cohorts/all_nr);
for ci = 1:n_cohorts
    ax(ci) = subplot(all_nr,all_nc,r2c(all_nr,all_nc,ci));
    hold on
    title(ax(ci),fgfr3_effects(ci))
end

for mut_ind = 1:2
    for ant_ind = 1:2
        t = arrayifyNonuniform(out_all{mut_ind,ant_ind},"t");
        max_t = max(t,[],'all');
        for i = 1:numel(out_all{mut_ind,ant_ind})
            out_all{mut_ind,ant_ind}(i).rint = sum(out_all{mut_ind,ant_ind}(i).avg.*(options.rr.^1)',1,'omitnan')./sum(out_all{mut_ind,ant_ind}(i).avg,1,'omitnan');
        end
        rint = arrayifyNonuniform(out_all{mut_ind,ant_ind},"rint");
        for ci = 1:n_cohorts
            [xx,yy,pc] = my_patchPlot(squeeze(t(ci,:,:))',squeeze(rint(ci,:,:))',true);
            color_ind = sub2ind([2,2],mut_ind,ant_ind);
            pp(ci,mut_ind,ant_ind) = patch(ax(ci),pc{1},max(0,pc{2}),type_colors(color_ind,:),'FaceAlpha',0.15,'EdgeColor','none','DisplayName',tumor_type(mut_ind,ant_ind));
            ll(ci,mut_ind,ant_ind) = plot(ax(ci),xx,yy,'Color',type_colors(color_ind,:),"LineWidth",1.25,'DisplayName',tumor_type(mut_ind,ant_ind));
            xlim(ax(ci),[0 max_t])
        end
    end
end
normalizeXLims(ax)
xlabel(ax,'Time (days)')
ylabel(ax,'r (cell widths)')
% title(["Radius of Maximal Clustering","of CTLs near HA Cells"])
set(ax,'FontSize',16)
L = legend(reshape(ll(1,:,:),[],1),"Location","best");
L.Title.String = "Tumor Type";
f(end).Position(3:4) = [703 420];
%% print
printFigures(f,cohort_name);

