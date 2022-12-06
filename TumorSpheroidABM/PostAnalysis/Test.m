clearvars;

%% load cohort file
cohort_files = dir("../data/cohort_*");
load(sprintf('%s/%s',cohort_files(end).folder,cohort_files(end).name))

%% load each sample file's tracked variable
for i = total_runs:-1:1
    load(sprintf("../data/%d/output_final.mat",ids(i)),"tracked")
    all_tracked(i) = tracked;
end

%% pull out t and NT
all_t = [all_tracked.t];
unique_t = unique(all_t);
all_NT = [all_tracked.NT];
unique_NT = unique(all_NT);

%% mesh
[tt,nn] = ndgrid(unique_t,unique_NT);
pp = zeros(size(tt));

%% determine largest tumor size prior to any time in unique_t
largest_previous_tumor_size = zeros([length(all_t),total_runs]);
for si = 1:total_runs
    for ti = 1:length(all_tracked(si).t)
        tinds = unique_t>=all_tracked(si).t(ti);
        largest_previous_tumor_size(tinds,si) = max(all_tracked(si).NT(ti),largest_previous_tumor_size(tinds,si));
    end
end

%% proportion of all samples that are size N by time T graphed at (T,N)
larger_than_n_by_t = false([size(tt),total_runs]);
for si = 1:total_runs
    for ti = 1:length(all_tracked(si).t)
        tinds = unique_t>=all_tracked(si).t(ti);
        ninds = unique_NT<=all_tracked(si).NT(ti);
        larger_than_n_by_t(tinds,ninds,si) = true;
    end
end

pp = mean(larger_than_n_by_t,3);

%%

figure;
mesh(tt,nn,pp)
xlabel('t (days)')
ylabel('Tumor (cells)')

%% variance for each n
tstar = zeros(total_runs,length(unique_NT));
for ni = 1:length(unique_NT)
    for si = 1:total_runs
        temp = unique_t(find(larger_than_n_by_t(:,ni,si),1));
        if isempty(temp)
            tstar(si,ni) = unique_t(end);
        else
            tstar(si,ni) = temp;
        end
    end
end

v = std(tstar,[],1);
v_normalized = std(tstar./mean(tstar,1),[],1);

[~,best_n_ind] = max(v);
best_n = unique_NT(best_n_ind);

%% time to best n
time_to_best_n = tstar(:,best_n_ind);
% for si = 1:total_runs
%     time_to_best_n(si) = all_tracked(si).t(find(all_tracked(si).NT>=best_n,1));
% end

%% get immune proximity probabilities
r = 2;
[P,T] = averageDistanceToImmuneCohort(sprintf('%s/%s',cohort_files(end).folder,cohort_files(end).name),r);


%%
figure;
hold on;

colors = flipud(cubehelix(255,3,-.1,1,1,[.2,.8],[.3,.7]));
c = colorbar;
colormap(colors)
caxis([min(time_to_best_n),max(time_to_best_n)])
for si = 1:total_runs
    color_ind = round(interp1([min(time_to_best_n),max(time_to_best_n)],[1,size(colors,1)],time_to_best_n(si)));
    plot(T(:,si),P(:,si),'Color',colors(color_ind,:),'LineWidth',1)
end
xlabel('Time (d)')
ylabel('Proportion')
title(["Proportion of Tumor Cells","With an Immune Cell Within r = 2"])
c.Label.String = sprintf("Time to %d Tumor Cells",best_n);
c.Label.Rotation = 270;
c.Label.VerticalAlignment = "baseline";
set(gca,'FontSize',20)
%% scatter the above proportions against time to size N
figure; hold on;
ylim([min(time_to_best_n),max(time_to_best_n)])
[pl(1),pl(2)] = bounds(P(:));
for ti = 1:size(T,1)
    if ti>1
        delete(sc_plot)
        delete(interp_line)
    end
    sc_plot = scatter(P(ti,:),time_to_best_n);

    beta = [P(ti,:)'.^2,P(ti,:)',ones(total_runs,1)]\time_to_best_n;

    xx = linspace(pl(1),pl(2),101);
    interp_line = plot(xx,beta(end-2)*xx.^2+beta(end-1)*xx + beta(end));

    title(sprintf("t = %3.2f days",T(ti,1)))
    drawnow
    pause(0.2)
end


