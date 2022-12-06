function f = plotEventByTumorType(tracked,fn,cohort_labels,nsamps_per_condition,name,options)

tumor_type = ["LA","HA";"LA + Mut","HA + Mut"];

%% get cohort info
n_cohorts = numel(cohort_labels);
f = gobjects(0,1);
all_nr = ceil(sqrt(n_cohorts));
all_nc = ceil(n_cohorts/all_nr);

t = arrayifyNonuniform(tracked,"t");
V = arrayifyNonuniform(tracked,fn);
colors = lines(n_cohorts);
type_colors = jet(4);

t = reshape(t,n_cohorts,nsamps_per_condition,[]);
V = reshape(V,n_cohorts,nsamps_per_condition,[],2,2);
tracked = reshape(tracked,[n_cohorts,nsamps_per_condition]);

if options.print_everything
    %% plot all events per day
    ax = gobjects([n_cohorts,2,2]);
    for mut_ind = 1:2
        for ant_ind = 1:2
            f(end+1,1) = figure;
            for ci = 1:n_cohorts

                ax(ci,mut_ind,ant_ind) = subplot(all_nr,all_nc,ci);
                ax(ci,mut_ind,ant_ind).XLim(2) = eps;
                hold on;
                for si = 1:nsamps_per_condition
                    plot(ax(ci,mut_ind,ant_ind),tracked(ci,si).t(1:end-1),tracked(ci,si).(fn)(2:end,mut_ind,ant_ind)./diff(tracked(ci,si).t));
                    ax(ci,mut_ind,ant_ind).XLim(2) = max(ax(ci,mut_ind,ant_ind).XLim(2),tracked(ci,si).t(end));
                end
                title(sprintf("%s with %s %s",tumor_type(mut_ind,ant_ind),cohort_labels(ci),name))


            end
            f(end).Position(3:4) = [560   257];
            normalizeYLims(f(end))
        end
    end

    xlabel(ax,'Time (days)')
    ylabel(ax,[sprintf("%s Events",name),"Per Day (#/d)"])

    %% patch plot events per day
    f(end+1,1) = figure;
    ax = gobjects(2,2);
    pp = gobjects(n_cohorts,2,2);
    ll = gobjects(n_cohorts,2,2);
    max_t = max(t,[],'all');
    for mut_ind = 1:2
        for ant_ind = 1:2
            ax(mut_ind,ant_ind) = subplot(2,2,r2c(2,2,[mut_ind,ant_ind]));
            hold on;
            for ci = 1:n_cohorts
                clear r
                for si = 1:nsamps_per_condition
                    r(si)=ksr(squeeze(t(ci,si,1:end-1)),squeeze(V(ci,si,2:end,mut_ind,ant_ind))./diff(squeeze(t(ci,si,:))),4/24,ceil(max_t*4));
                end
                [xx,yy,pc] = my_patchPlot(arrayify(r,"x",1),arrayify(r,"f",1),true);
                pp(ci,mut_ind,ant_ind) = patch(pc{1},max(0,pc{2}),colors(ci,:),'FaceAlpha',0.15,'EdgeColor','black');
                ll(ci,mut_ind,ant_ind) = plot(xx,yy,'Color',colors(ci,:),"LineWidth",1.25);
            end
            title(ax(mut_ind,ant_ind),sprintf("%s %s",tumor_type(mut_ind,ant_ind),name))
        end
    end

    f(end).Position(3:4) = [1169 388];
    L = legend(ax(1,1),reshape(ll(:,1,1),[],1),cohort_labels(:),"Location","northwest");
    set(ax,'FontSize',16)
    xlabel(ax,'Time (days)')
    ylabel(ax,[sprintf("%s Events",name),"Per Day (#/d)"])

    %% all rates
    ax = gobjects(n_cohorts,2,2);
    for mut_ind = 1:2
        for ant_ind = 1:2
            f(end+1,1) = figure;
            for ci = 1:n_cohorts
                ax(ci,mut_ind,ant_ind) = subplot(all_nr,all_nc,ci);
                ax(ci,mut_ind,ant_ind).XLim(2) = eps;
                hold on;
                for si = 1:nsamps_per_condition
                    plot(ax(ci,mut_ind,ant_ind),tracked(ci,si).t(1:end-1),tracked(ci,si).(fn)(2:end,mut_ind,ant_ind)./(diff(tracked(ci,si).t).*tracked(ci,si).tumor_types(1:end-1,mut_ind,ant_ind)))
                    ax(ci,mut_ind,ant_ind).XLim(2) = max(ax(ci,mut_ind,ant_ind).XLim(2),tracked(ci,si).t(end));
                end
                title(sprintf("%s with %s %s",tumor_type(mut_ind,ant_ind),cohort_labels(ci),name))

            end
            normalizeYLims(f(end))
            f(end).Position(3:4) = [560   257];
        end
    end

    xlabel(ax,'Time (days)')
    ylabel(ax,sprintf("%s Rate (d^{-1})",name))
end

TumTypes = arrayifyNonuniform(tracked,"tumor_types");
V_per = V(:,:,2:end,:,:)./TumTypes(:,:,1:end-1,:,:);

%% patch rates
f(end+1,1) = figure("Name",sprintf("%s_rates",fn));
ax = gobjects(2,2);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);
max_t = max(t,[],'all');
for mut_ind = 1:2
    for ant_ind = 1:2
        ax(mut_ind,ant_ind) = subplot(2,2,r2c(2,2,[mut_ind,ant_ind]));
        hold on;
        for ci = 1:n_cohorts
            clear r
            for si = 1:nsamps_per_condition
                r(si)=ksr(squeeze(t(ci,si,1:end-1)),squeeze(V_per(ci,si,:,mut_ind,ant_ind))./diff(squeeze(t(ci,si,:))),4/24,ceil(max_t*4));
            end
            [xx,yy,pc] = my_patchPlot(arrayify(r,"x",1),arrayify(r,"f",1),true);
            pp(ci,mut_ind,ant_ind) = patch(pc{1},max(0,pc{2}),colors(ci,:),'FaceAlpha',0.15,'EdgeColor','black');
            ll(ci,mut_ind,ant_ind) = plot(xx,yy,'Color',colors(ci,:),"LineWidth",1.25);
        end
        title(ax(mut_ind,ant_ind),sprintf("%s %s",tumor_type(mut_ind,ant_ind),name))
    end
end

normalizeXLims(f(end))
normalizeYLims(f(end))
f(end).Position(3:4) = [1263 516];
L = legend(ax(1,1),reshape(ll(:,1,1),[],1),cohort_labels(:),"Location","northwest");
set(ax,'FontSize',16)
xlabel(ax,'Time (days)')
ylabel(ax,sprintf("%s Rate (d^{-1})",name))

%% patch rates by condition
f(end+1,1) = figure("Name",sprintf("%s_rates_by_condition",fn));
ax = gobjects(n_cohorts,1);
pp = gobjects(n_cohorts,2,2);
ll = gobjects(n_cohorts,2,2);
max_t = max(t,[],'all');
for ci = 1:n_cohorts
    ax(ci) = subplot(all_nr,all_nc,r2c(all_nr,all_nc,ci));
    hold on;
    for mut_ind = 1:2
        for ant_ind = 1:2
            clear r
            for si = 1:nsamps_per_condition
                r(si)=ksr(squeeze(t(ci,si,1:end-1)),squeeze(V_per(ci,si,:,mut_ind,ant_ind))./diff(squeeze(t(ci,si,:))),4/24,ceil(max_t*4));
            end
            [xx,yy,pc] = my_patchPlot(arrayify(r,"x",1),arrayify(r,"f",1),true);
            pp(ci,mut_ind,ant_ind) = patch(pc{1},max(0,pc{2}),type_colors(mut_ind+(ant_ind-1)*2,:),'FaceAlpha',0.15,'EdgeColor','black');
            ll(ci,mut_ind,ant_ind) = plot(xx,yy,'Color',type_colors(mut_ind+(ant_ind-1)*2,:),"LineWidth",1.25);
        end
    end
    title(ax(ci),sprintf("%s %s",cohort_labels(ci),name))
end

normalizeXLims(f(end))
normalizeYLims(f(end))
f(end).Position(3:4) = [1263 516];
L = legend(ax(1),reshape(ll(1,:,:),[],1),tumor_type(:),"Location","northwest");
set(ax,'FontSize',16)
xlabel(ax,'Time (days)')
ylabel(ax,sprintf("%s Rate (d^{-1})",name))
