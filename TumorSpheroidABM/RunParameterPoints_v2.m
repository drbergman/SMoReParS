% this script will run samples across many parameter vectors where the
% parameter vectors are in a lattice_parameters type struct array selected
% by smorepars. it will also look at parameter vectors outside of this set
% to see if the rejected parameters produce a poor fit with the data

% this script actually does not currently run samples, but rather uses the
% original 3^7 grid and determines which of these are accepted based on
% smorepars, looking at those trajectories for each of their 6 samples. It
% looks at the RSS for these grouped by acceptance/rejection

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

%% cohort structure
cohort_pars.nsamps_per_condition = 6;
cohort_pars.min_parfor_num = 4;

%%
M = allBaseParameters();
%%

M.setup.ndims = 2;
M.setup.censor_date = 3;
M.setup.N0 = 1e2;
M.setup.agent_initialization_location = "uniform";
M.setup.use_rates_for_intitial_proportions = false;

M.save_pars.make_save = false;
M.save_pars.dt = Inf;

M.pars.max_dt = 0.25 / 24; % number of days per step
M.pars.apop_rate = 0;

M.cycle_pars.dna_check_g1 = false;
M.cycle_pars.dna_check_s = false;
M.cycle_pars.dna_check_g2 = false;
M.cycle_pars.dna_check_m = false;

M.cycle_pars.arrest_prob_g1 = 0.05;
M.cycle_pars.arrest_prob_s = 0.00;
M.cycle_pars.arrest_prob_g2 = 0.05;
M.cycle_pars.arrest_prob_m = 0.00;

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = false;

%% load the ABM parameters using the best fit values of lambda and alpha; simulate

load("ProfileLikelihood/ABMParamEstimates_FromProfile_WithK.mat","LP1","abm_region_1_log")
% for i = 1:numel(LP1(1).values)
%     for vpi = 1:numel(LP1)
%         M = setField(M,LP1(vpi).path,LP1(vpi).values(i,:));
%     end
%     for j = 1:cohort_pars.nsamps_per_condition
%         out(j,i) = simPatient(M);
%     end
% end

%% load the ABM parameters using the above and intersecting with the corresponding region from the best fit value of K
% load("ProfileLikelihood/ABMParamEstimates.mat","LP2")
% also_after_restriction = ismember([LP1.values],[LP2.values],"rows");

%% compare ABM output with data

% Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
% Data point for 5 hours taken out, since it is incommensurate.

tt = [0      10     24     36     48     72    ]';        % hours
tt = tt/24;                                             % days

% Control data
data = [0.899  1.340  1.633  2.408  3.557  5.583]';   % millions of cells
data_std = [0.099  0.193  0.207  0.298  0.168  0.364]';   % millions of cells

factor = 100/data(1);

data = data*factor;
data_std = data_std*factor;

figure;
ax = gobjects(2,1);
for ri = 1:2
    ax(ri) = subplot(2,1,r2c(2,1,[ri,1])); hold on;
end
patch(ax(1),[tt;flip(tt)],[data-data_std;flip(data+data_std)],"blue","FaceAlpha",0.2,"EdgeColor","none")
exp_plot = plot(ax(1),tt,data,"blue","LineStyle","--","LineWidth",2,"Marker","o","DisplayName","Exp Data");
title(ax(1),"Restricted by \lambda and \alpha")

set(ax,"FontSize",20)
ylabel(ax(1),"Cell Count")
xlabel(ax(1),"Time (d)")

%% plot results that do match these parameters 
C = load("data/cohort_230124175743017/output.mat","ids","lattice_parameters","cohort_size");
sz = size(C.ids);
vpi = cell(7,1);
par_vals = zeros(1,7);
y_temp = [];
use_abm_region_1_var = isequal(C.cohort_size,size(abm_region_1_log));
for i = 1:numel(C.ids)
    [vpi{:},si] = ind2sub(sz,i);
    if use_abm_region_1_var
        accepted = abm_region_1_log(vpi{:});
    else
        for j = 1:7
            par_vals(j) = C.lattice_parameters(j).values(vpi{j});
        end
        accepted = ismember(par_vals,[LP1.values],"rows"); % If I just saved the logical array that I used to pick the abm pars, then use that
    end
    if accepted
        load(sprintf("data/sims/%s/output_final.mat",C.ids(vpi{:},si)),"tracked")
        y_temp = [y_temp,tracked.NT];
%         abm_line(i) = plot(ax(1),tracked.t,tracked.NT,"Color",[0 0 0 0.1],"DisplayName","ABM Sample");
    end
end
[x,y_mean,pc] = my_patchPlot(tracked.t,y_temp);
patch(ax(1),pc{1},pc{2},"green","FaceAlpha",0.2,"EdgeColor","none");
abm_mean = plot(ax(1),x,y_mean,"green","LineWidth",2,"DisplayName","ABM Mean");

% legend(ax(1),[exp_plot(1);abm_mean;abm_line(find(isgraphics(abm_line),1))],"Location","best")
legend(ax(1),[exp_plot(1);abm_mean],"Location","best")

%% histograms out objfn for those selected/those not
sz = size(C.ids);
vpi = cell(7,1);
par_vals = zeros(1,7);
H = cell(1,2);
tt_min = round(1440*tt);
% F = @(t) sum(((interp1(round(t.t*1440),t.NT,tt_min)-data)./data_std).^2,'all');
F = @(t) sum(((interp1(round(t.t*1440),t.NT,tt_min)-data)).^2,'all');
use_abm_region_1_var = isequal(C.cohort_size,size(abm_region_1_log));

for i = 1:numel(C.ids)
    [vpi{:},si] = ind2sub(sz,i);
    load(sprintf("data/sims/%s/output_final.mat",C.ids(i)),"tracked")
    if use_abm_region_1_var
        accepted = abm_region_1_log(vpi{:});
    else
        for j = 1:7
            par_vals(j) = C.lattice_parameters(j).values(vpi{j});
        end
        accepted = ismember(par_vals,[LP1.values],"rows"); % If I just saved the logical array that I used to pick the abm pars, then use that
    end
    val = F(tracked);
    if ~accepted
        H{2} = [H{2};val];
    else
        H{1} = [H{1};val];
    end
end

%% actually plot the above histograms
histogram(ax(2),H{1},"FaceColor","green","Normalization","count","DisplayName","Selected Samples")
histogram(ax(2),H{2},"FaceColor","red","Normalization","count","DisplayName","Rejected Samples")

legend(ax(2))
xlabel(ax(2),"RSS")
ylabel(ax(2),"Count")

%%
% savefig("ItWorked")

%% residuals at time points
sz = size(C.ids);
vpi = cell(7,1);
par_vals = zeros(1,7);
R = cell(1,2);
tt_min = round(1440*tt);
% F = @(t) sum(((interp1(round(t.t*1440),t.NT,tt_min)-data)./data_std).^2,'all');
F = @(t) sum(((interp1(round(t.t*1440),t.NT,tt_min)-data)).^2,'all');
use_abm_region_1_var = isequal(C.cohort_size,size(abm_region_1_log));

for i = 1:numel(C.ids)
    load(sprintf("data/sims/%s/output_final.mat",C.ids(i)),"tracked")
    [vpi{:},si] = ind2sub(sz,i);
    if use_abm_region_1_var
        accepted = abm_region_1_log(vpi{:});
    else
        for j = 1:7
            par_vals(j) = C.lattice_parameters(j).values(vpi{j});
        end
        accepted = ismember(par_vals,[LP1.values],"rows"); % If I just saved the logical array that I used to pick the abm pars, then use that
    end
    res_temp = interp1(round(tracked.t*1440),tracked.NT,tt_min)-data;
    if ~accepted
        R{2} = [R{2},res_temp];
    else
        R{1} = [R{1},res_temp];
    end
end

%% plot residuals from above
figure;
ax = gobjects(length(tt),1);
for j = 1:length(tt)
    min_temp = min(min(R{1}(j,:)),min(R{2}(j,:)));
    max_temp = max(max(R{1}(j,:)),max(R{2}(j,:)));

    ax(j) = subplot(1,length(tt),j); hold on;
    [~,binEdges] = histcounts([R{1}(j,:),R{2}(j,:)]./data(j));
    histogram(ax(j),R{1}(j,:)./data(j),"BinEdges",binEdges,"FaceColor","green","Normalization","pdf","DisplayName","Selected Samples","EdgeColor","none")
    histogram(ax(j),R{2}(j,:)./data(j),"BinEdges",binEdges,"FaceColor","red","Normalization","pdf","DisplayName","Rejected Samples","EdgeColor","none")
    if min_temp~=0 || max_temp~=0
        xlim(max(abs([min_temp,max_temp])) * [-1 1]/data(j))
    end
    title(sprintf("t = %3.2f d",tt(j)))
    xlabel("Res / data")
end
legend(ax(2))

%%
rmpath("~/Documents/MATLAB/myfunctions/")