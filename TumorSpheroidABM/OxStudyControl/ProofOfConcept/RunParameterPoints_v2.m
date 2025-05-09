% this script will run samples across many parameter vectors where the
% parameter vectors are in a lattice_parameters type struct array admitted
% by smorepars. it will also look at parameter vectors outside of this set
% to see if the rejected parameters produce a poor fit with the data

% this script actually does not currently run samples, but rather uses the
% original 3^7 grid and determines which of these are admitted based on
% smorepars, looking at those trajectories for each of their 6 samples. It
% looks at the RSS for these grouped by admittance/rejection

clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions

%% load the ABM parameters using the best fit values of lambda, alpha, and K; identify the admitted/accepted pars
load("../../OxStudyControl/ProfileLikelihood/data/ABMParamEstimates_FromProfile_WithK.mat","LP1","abm_region_1_log")
C = load("../../data/cohort_230124175743017/output.mat","ids","lattice_parameters","cohort_size","nsamps_per_condition");
sz = size(C.ids);
vpi = cell(7,1);
par_vals = zeros(1,7);
use_abm_region_1_var = isequal(C.cohort_size,size(abm_region_1_log));
if use_abm_region_1_var
    admitted = abm_region_1_log(:);
else
    admitted = false(numel(C.ids),1);
    for i = 1:numel(C.ids)
        [vpi{:},si] = ind2sub(sz,i);
        for j = 1:7
            par_vals(j) = C.lattice_parameters(j).values(vpi{j});
        end
        admitted(i) = ismember(par_vals,[LP1.values],"rows"); % If I just saved the logical array that I used to pick the abm pars, then use that
    end
end
admitted_ind = find(admitted);
rejected_ind = find(~admitted);

%% load and plot experimental data
load("../../OxStudyControl/ODEFitting/data/ExperimentalData.mat")

figure;
ax = gobjects(2,1);
for ri = 1:2
    ax(ri) = subplot(2,1,r2c(2,1,[ri,1])); hold on;
end
patch(ax(1),[tt;flip(tt)],[count-count_std;flip(count+count_std)],"blue","FaceAlpha",0.2,"EdgeColor","none")
exp_plot = plot(ax(1),tt,count,"blue","LineStyle","--","LineWidth",2,"Marker","o","DisplayName","Exp Data");
title(ax(1),"Restricted by \lambda, \alpha, and K")

set(ax,"FontSize",20)
ylabel(ax(1),"Cell Count")
xlabel(ax(1),"Time (d)")

%% plot results that do match these parameters 
load(sprintf("../../data/sims/%s/output_final.mat",C.ids(1)),"tracked")
nt = length(tracked.t);
y_temp = zeros(nt,length(admitted_ind));
C.ids = reshape(C.ids,[],C.nsamps_per_condition);
for i = 1:length(admitted_ind)
    for j = 1:C.nsamps_per_condition
        load(sprintf("../../data/sims/%s/output_final.mat",C.ids(admitted_ind(i),j)),"tracked")
        y_temp(:,i,j) = tracked.NT;
    end
end
y_temp = reshape(y_temp,nt,[]);
[x,y_mean,pc] = patchPlotCoords(tracked.t,y_temp);
patch(ax(1),pc{1},pc{2},"green","FaceAlpha",0.2,"EdgeColor","none");
abm_mean = plot(ax(1),x,y_mean,"green","LineWidth",2,"DisplayName","ABM Mean");

% legend(ax(1),[exp_plot(1);abm_mean;abm_line(find(isgraphics(abm_line),1))],"Location","best")
legend(ax(1),[exp_plot(1);abm_mean],"Location","best")

%% histograms out objfn for those admitted/those not
H = {zeros(length(admitted_ind),C.nsamps_per_condition),zeros(length(rejected_ind),C.nsamps_per_condition)};
tt_min = round(1440*tt);
% F = @(t) sum(((interp1(round(t.t*1440),t.NT,tt_min)-count)).^2,'all'); % round saved data time points to nearest minute (to prevent, for example, final timepoint being a little smaller than 3)
F = @(t) sum(((interp1(round(t.t*1440),t.NT,tt_min)-count)./count_std).^2,'all'); % round saved data time points to nearest minute (to prevent, for example, final timepoint being a little smaller than 3)


for i = 1:length(admitted_ind)
    for j = 1:C.nsamps_per_condition
        load(sprintf("../../data/sims/%s/output_final.mat",C.ids(admitted_ind(i),j)),"tracked")
        H{1}(i,j) = F(tracked);
    end
end
H{1} = H{1}(:)';
for i = 1:length(rejected_ind)
    for j = 1:C.nsamps_per_condition
        load(sprintf("../../data/sims/%s/output_final.mat",C.ids(rejected_ind(i),j)),"tracked")
        H{2}(i,j) = F(tracked);
    end
end
H{2} = H{2}(:)';

%% actually plot the above histograms
[~,binEdges] = histcounts([H{:}]);
histogram(ax(2),H{1},binEdges,"FaceColor","green","Normalization","count","DisplayName","Admitted Samples")
histogram(ax(2),H{2},binEdges,"FaceColor","red","Normalization","count","DisplayName","Rejected Samples")

legend(ax(2))
xlabel(ax(2),"nRSS")
ylabel(ax(2),"Count")

savefig("figures/fig/AdmittedSamples.fig")
print("figures/png/AdmittedSamples.png","-dpng")

%%
% savefig("ItWorked")

%% residuals at time points
R = {zeros(length(tt),length(admitted_ind),C.nsamps_per_condition),zeros(length(tt),length(rejected_ind),C.nsamps_per_condition)};
tt_min = round(1440*tt);
for i = 1:length(admitted_ind)
    for j = 1:C.nsamps_per_condition
        load(sprintf("../../data/sims/%s/output_final.mat",C.ids(admitted_ind(i),j)),"tracked")
        R{1}(:,i,j) = interp1(round(tracked.t*1440),tracked.NT,tt_min)-count;
    end
end
R{1} = reshape(R{1},length(tt),[]);
for i = 1:length(rejected_ind)
    for j = 1:C.nsamps_per_condition
        load(sprintf("../../data/sims/%s/output_final.mat",C.ids(rejected_ind(i),j)),"tracked")
        R{2}(:,i,j) = interp1(round(tracked.t*1440),tracked.NT,tt_min)-count;
    end
end
R{2} = reshape(R{2},length(tt),[]);

%% plot residuals from above
figure;
ax = gobjects(length(tt),1);
for j = 1:length(tt)
    min_temp = min(min(R{1}(j,:)),min(R{2}(j,:))); % minimum residual at time tt(j)
    max_temp = max(max(R{1}(j,:)),max(R{2}(j,:))); % maximum residual at time tt(j)

    ax(j) = subplot(1,length(tt),j); hold on;
    [~,binEdges] = histcounts([R{1}(j,:),R{2}(j,:)]./count_std(j));
    histogram(ax(j),R{1}(j,:)./count_std(j),"BinEdges",binEdges,"FaceColor","green","Normalization","pdf","DisplayName","Admitted Samples","EdgeColor","none")
    histogram(ax(j),R{2}(j,:)./count_std(j),"BinEdges",binEdges,"FaceColor","red","Normalization","pdf","DisplayName","Rejected Samples","EdgeColor","none")
    if min_temp~=0 || max_temp~=0
        xlim(max(abs([min_temp,max_temp])) * [-1 1]/count_std(j)) % set xlim to be big enough to get largest |residual| and center at x=0
    end
    title(sprintf("t = %3.2f d",tt(j)))
    xlabel("Z Score")
end
legend(ax(2))

savefig("figures/fig/TimeSeriesZScores_ByAdmittance.fig")
print("figures/png/TimeSeriesZScores_ByAdmittance.png","-dpng")