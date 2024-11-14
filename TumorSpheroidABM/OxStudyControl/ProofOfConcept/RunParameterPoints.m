% this script will run samples across many parameter vectors where the
% parameter vectors are in a lattice_parameters type struct array selected
% by smorepars. it will also look at parameter vectors outside of this set
% to see if the rejected parameters produce a poor fit with the data

clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../..")

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

M.chemo_pars.dna_check_g1 = false;
M.chemo_pars.dna_check_s = false;
M.chemo_pars.dna_check_g2 = false;
M.chemo_pars.dna_check_m = false;

M.chemo_pars.arrest_coeff_g1 = 0.05;
M.chemo_pars.arrest_coeff_s = 0.00;
M.chemo_pars.arrest_coeff_g2 = 0.05;
M.chemo_pars.arrest_coeff_m = 0.00;

M.plot_pars.plot_fig = false;
M.plot_pars.plot_location = false;

%% load the ABM parameters using the best fit values of lambda and alpha; simulate

load("../ProfileLikelihood/data/ABMParamEstimates.mat","LP1")
for i = 1:numel(LP1(1).values)
    for vpi = 1:numel(LP1)
        M = setField(M,LP1(vpi).path,LP1(vpi).values(i,:));
    end
    for j = 1:cohort_pars.nsamps_per_condition
        out(j,i) = simPatient(M);
    end
end

%% load the ABM parameters using the above and intersecting with the corresponding region from the best fit value of K
load("../ProfileLikelihood/data/ABMParamEstimates.mat","LP2")
also_after_restriction = ismember([LP1.values],[LP2.values],"rows");

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
ax = gobjects(2,2);
for ri = 1:2
    for ci = 1:2
    ax(ri,ci) = subplot(2,2,r2c(2,2,[ri,ci])); hold on;
    end
end
for i = 1:numel(out)
    abm_line(i) = plot(ax(1,1),out(i).tracked.t,out(i).tracked.NT,"Color","black","DisplayName","ABM Sample");
    [ri,~] = ind2sub(size(out),i);
    if also_after_restriction(ri) % for simulations that also came from using the fitted value of K
        plot(ax(1,2),out(i).tracked.t,out(i).tracked.NT,"Color","black")
    end
end
for axi = 1:2
    patch(ax(1,axi),[tt;flip(tt)],[data-data_std;flip(data+data_std)],"blue","FaceAlpha",0.2,"EdgeColor","none")
    exp_plot(axi) = plot(ax(1,axi),tt,data,"blue","LineStyle","--","LineWidth",2,"Marker","o","DisplayName","Exp Data");
end
title(ax(1,1),"Restricted by \lambda and \alpha")
title(ax(1,2),"Restricted ALSO by K")

%% extract sim data
NT = zeros([size(out(1).tracked.t,1),size(out)]);
for i = 1:size(out,1)
    for j = 1:size(out,2)
        NT(:,i,j) = out(i,j).tracked.NT;
    end
end
%% plot the abm mean +/- SD for the two restrictions chosen
[x,y_mean,pc] = patchPlotCoords(out(1).tracked.t,reshape(NT,size(NT,1),[]));
patch(ax(1,1),pc{1},pc{2},"green","FaceAlpha",0.2,"EdgeColor","none")
abm_plot(1) = plot(ax(1,1),x,y_mean,"green","LineWidth",2,"DisplayName","ABM mean");

[x,y_mean,pc] = patchPlotCoords(out(1).tracked.t,reshape(NT(:,:,also_after_restriction),size(NT,1),[]));
patch(ax(1,2),pc{1},pc{2},"green","FaceAlpha",0.2,"EdgeColor","none")
abm_plot(1,2) = plot(ax(1,2),x,y_mean,"green","LineWidth",2,"DisplayName","ABM mean");

set(ax,"FontSize",20)
legend(ax(1,1),[exp_plot(1);abm_plot(1);abm_line(1)],"Location","best")
ylabel(ax(1,1),"Cell Count")
xlabel(ax,"Time (d)")

%% plot results that do no match these parameters (this was unhelpful for determining if the rejected parameters fit the data)
% C = load("data/cohort_230124175743017/output.mat","ids","lattice_parameters");
% sz = size(C.ids);
% vpi = cell(7,1);
% par_vals = zeros(1,7);
% for i = 1:numel(C.ids)
%     [vpi{:},si] = ind2sub(sz,i);
%     for j = 1:7
%         par_vals(j) = C.lattice_parameters(j).values(vpi{j});
%     end
%     if ~ismember(par_vals,[LP1.values],"rows")
%         load(sprintf("data/sims/%s/output_final.mat",C.ids(vpi{:},si)),"tracked")
%         plot(ax(2,1),tracked.t,tracked.NT,"Color","black")
%     end
%     if ~ismember(par_vals,[LP2.values],"rows")
%         load(sprintf("data/sims/%s/output_final.mat",C.ids(vpi{:},si)),"tracked")
%         plot(ax(2,2),tracked.t,tracked.NT,"Color","black")
%     end
% end

%% histograms out objfn for those selected/those not
C = load("../../data/cohort_230124175743017/output.mat","ids","lattice_parameters");
sz = size(C.ids);
vpi = cell(7,1);
par_vals = zeros(1,7);
H = cell(2,2);
tt_min = round(1440*tt);
% F = @(t) sum(((interp1(round(t.t*1440),t.NT,tt_min)-data)./data_std).^2,'all');
F = @(t) sum(((interp1(round(t.t*1440),t.NT,tt_min)-data)).^2,'all');
for i = 1:numel(C.ids)
    [vpi{:},si] = ind2sub(sz,i);
    load(sprintf("../../data/sims/%s/output_final.mat",C.ids(i)),"tracked")
    for j = 1:7
        par_vals(j) = C.lattice_parameters(j).values(vpi{j});
    end
    val = F(tracked);
    if ~ismember(par_vals,[LP1.values],"rows")
        H{1,2} = [H{1,2};val];
    else
        H{1,1} = [H{1,1};val];
    end
    if ~ismember(par_vals,[LP2.values],"rows")
        H{2,2} = [H{2,2};val];
    else
        H{2,1} = [H{2,1};val];
    end
end

%% actually plot the above histograms
for i = 1:2 % the two restrictions
    histogram(ax(2,i),H{i,1},"FaceColor","green","Normalization","pdf","DisplayName","Selected Samples")
    histogram(ax(2,i),H{i,2},"FaceColor","red","Normalization","pdf","DisplayName","Rejected Samples")
end
legend(ax(2,1))
xlabel(ax(2,:),"RSS")
ylabel(ax(2,:),"PDF")

%% save figure
% savefig("ItWorked")

%%
rmpath("../..")

