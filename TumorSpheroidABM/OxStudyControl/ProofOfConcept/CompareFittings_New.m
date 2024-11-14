% This script will compare fitting the ABM using SMoRe ParS against using
% the data directly.

% This is a version of CompareFittings.m that is being updated to use the
% new workflow.

clearvars;
addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
cohort_name = "cohort_230124175743017";
files.abm_profile = "../SampleABM/data/AcceptedParameters_New_all_profiles_resampled.mat";
files.experimental_data = "../ODEFitting/data/ExperimentalData_New.mat";

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.resolution = '-r1200';

%% colors
accept_color = [0    0.4470    0.7410];
reject_color = [0.8500    0.3250    0.0980];
all_color = [0 0 0];
most_likely_color = [255 203 5]/255;
least_likely_color = [0 39 76]/255;
%%
files.abm_cohort = sprintf("../../data/%s/output.mat",cohort_name);

%% load SMoRe ParS-derived fitting
load(files.abm_profile,"accepted_parameters")
abm_data_file = load(files.abm_profile,"files");
abm_data_file = abm_data_file.files.abm_data;

%% load cohort data and experimental data
C = load(files.abm_cohort,"ids","cohort_size","nsamps_per_condition");
E = load(files.experimental_data);
nt = length(E.t);

%% load residuals of mean
load("data/ResidualsOfMean_New.mat")

%% Compute log-likelihood values
LL = -0.5*(nt*log(2*pi()) + sum(E.D.S.^2) + sum((residuals_of_mean./E.D.S).^2,1));

%% stratify LL by SMoRe ParS Acceptance/Rejection
LL_accepted = LL(:,accepted_parameters);
LL_rejected = LL(:,~accepted_parameters);
LL = LL(:);
%% plot histograms of both
figureOnRight; hold on;
ax = gca;
[~,binEdges] = histcounts(LL);
histogram(LL_accepted,"BinEdges",binEdges,"FaceColor",accept_color,"Normalization","pdf","DisplayName","Accepted Samples","EdgeColor","none")
histogram(LL_rejected,"BinEdges",binEdges,"FaceColor",reject_color,"Normalization","pdf","DisplayName","Rejected Samples","EdgeColor","none")
legend(ax,"Location","northwest","FontSize",16)
xlabel("log-likelihood","FontSize",16)
ylabel("PDF","FontSize",16)

%% distribution relative likelihood to best setup
N = sum(accepted_parameters,"all");
divider = 1e-2;
minRelProb = exp(min(LL)-max(LL));
lbase = divider/minRelProb;
fn = @logLinScaleX;
LL_sort = sort(LL,"ascend");
XData = {exp(LL_sort-max(LL)),exp(sort(LL_accepted)-max(LL)),exp(sort(LL_rejected)-max(LL))};
CDF = {linspace(0,1,numel(LL)),linspace(0,1,length(LL_accepted)),linspace(0,1,length(LL_rejected))};
Names = ["All","Accepted","Rejected"];
Colors = containers.Map(1:3,{all_color,accept_color,reject_color});

%% plot cdf
f = figureOnRight("Name","RelativeLikelihoodCDF_New");
ax = gca;
hold on;
for i = 1:3
    plot(fn(XData{i},divider,lbase),CDF{i},"DisplayName",Names(i),"Color",Colors(i),"LineWidth",0.5)
end
plot(fn(exp(LL_sort(1:end-N)-LL_sort(end)),divider,lbase),linspace(0,1,length(LL_rejected)),"DisplayName","Least Likely","Color",0.5*[1 1 1],"LineWidth",0.5,"LineStyle",":")
plot(fn(exp(LL_sort(end-N+1:end)-LL_sort(end)),divider,lbase),linspace(0,1,length(LL_accepted)),"DisplayName","Most Likely","Color",0*[1 1 1],"LineWidth",0.5,"LineStyle",":")
xlabel("Relative Likelihood","FontSize",16)
ylabel("CDF","FontSize",16)
% legend(gca,"location","west","FontSize",16,"AutoUpdate","off")
xL = xlim;
xL(1) = fn(minRelProb,divider,lbase);
minPow10 = log10(minRelProb);
maxPow10 = log10(divider);
lowerPow10 = flip(maxPow10:ceil(0.5*(minPow10-maxPow10)):minPow10);
lowerXTicks = fn(10.^lowerPow10,divider,lbase);
upperXTicks = .5:.5:1;
xticks(unique([lowerXTicks,upperXTicks]))
% xticklabels([divider*lbase.^(lowerXTicks-divider),upperXTicks])
xticklabels(1:numel(xticks))
xL(2) = 1;
xlim(xL)
xline(divider,"LineStyle","-.","LineWidth",0.5)
% annotation("textarrow",[.2875,.1875],[0.9,0.9],"String","Log Scale","FontSize",6);
% annotation("arrow",[0.3675,0.4675],[0.9,0.9]);
% annotation("textarrow",[.7625,.8625],[0.15,0.15],"String","Linear Scale","FontSize",6);
% annotation("arrow",[0.6675,0.5675],[0.15,0.15]);

f.Units = "inches";
f.Position(3) = 5/3;
f.Position(4) = 1;

set(ax,"FontSize",8)
% xtl = {'10^{-50}','10^{-26}','0.01','0.5','1'};
% xticklabels(xtl)
% ax.XAxis.TickLabelRotation = 0;

margin = struct("left",[],"right",0.02,"top",.05,"bottom",.29);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

saveFigures(f,save_fig_opts)

%% attempt to plot PDF of above
figureOnRight;
hold on;
for i = 1:3
    y = diff(CDF{i})./diff(XData{i});
    plot(fn(.5*(XData{i}(1:end-1)+XData{i}(2:end)),divider,lbase),y/max(y),"DisplayName",Names(i),"Color",Colors(i),"LineWidth",2)
end
set(gca,"YScale","log")
xlabel("Relative Likelihood","FontSize",16)
ylabel("PDF","FontSize",16)
legend(gca,"location","best","FontSize",16,"AutoUpdate","off")

%% attempt #2 to plot PDF of above
f = figureOnRight("Name","RelativeLikelihoodPDF_New");
hold on;
% x = [0,logspace(-30,0,20)];
x = linspace(0,1,26);
for i = 1:3
    temp = histcounts(XData{i},x,"Normalization","pdf");
    plot(.5*(x(1:end-1)+x(2:end)),temp,"DisplayName",Names(i),"Color",Colors(i),"LineWidth",2)
end
set(gca,"YScale","log")
xlabel("Relative Likelihood","FontSize",16)
ylabel("PDF","FontSize",16)
legend(gca,"location","best","FontSize",16,"AutoUpdate","off")

saveFigures(f,save_fig_opts)

%% plotting percentage of N best LLs that were accepted
f = figureOnRight("Name","AcceptanceOfMostLikely_New");
[~,order] = sort(LL,"descend");
accepted_sort = accepted_parameters(order);
y = zeros(2);
y(1,1) = sum(accepted_sort(1:N))/N; % percentage of accepted that are in the most likely category
y(2,1) = sum(~accepted_sort(1:N))/(numel(LL)-N); % percentage of rejected that are in the most likely category
y(1,2) = sum(accepted_sort(N+1:end))/N; % percentage of accepted that are in the most likely category
y(2,2) = sum(~accepted_sort(N+1:end))/(numel(LL)-N); % percentage of rejected that are in the most likely category
% y(:,1) = [sum(accepted_sort(1:N)),sum(~accepted_sort(1:N))]/N;
% y(:,2) = [sum(accepted_sort(N+1:end)),sum(~accepted_sort(N+1:end))]/(numel(LL)-N);
bar(categorical(["Accepted","Rejected"],["Accepted","Rejected"]),100*y,"stacked")
ylabel("Percentage")
set(gca,"FontSize",16)
title("SMoRe ParS Acceptance of Most Likely Parameters")
legend(["Most Likely","Least Likely"],"Location","northwest")

saveFigures(f,save_fig_opts)

%% reorient above to have categories be most likely and least likely and compare 
f = figureOnRight("Name","AcceptanceByLikelihoodQuantile_New");
ax=gca;
[LL_sort,order] = sort(LL,"ascend");
accepted_sort = accepted_parameters(order);
n_quants = 17;
% Q = quantile(LL_sort,[1-N/numel(LL),1]);
Q = quantile(LL_sort,n_quants);
I_prev = 1;
y = zeros(n_quants,2);
for i = 1:n_quants
    I = find(LL_sort<=Q(i),1,"last");
    y(i,1) = sum(accepted_sort(I_prev:I))/(I-I_prev+1);
    y(i,2) = 1-y(i,1);
    I_prev = I;
end
b=bar(linspace(0,100,n_quants),100*y,"stacked","BarWidth",1,"EdgeColor","none");
b(1).FaceColor = accept_color;
b(2).FaceColor = reject_color;
xlim([0,100]+50/(n_quants-1)*[-1 1])
ylim([0 100])
% legend(["Accepted","Rejected"],"location","northwest","AutoUpdate","off")
ylabel("Probability")
xlabel("Likelihood Percentile")
set(ax,"FontSize",16)
yticks(ax,0:25:100)
ax.XAxis.TickLength = [0 0];
xline(100*(1-N/numel(LL)),"LineWidth",0.5)

f.Units = "inches";
f.Position(3) = 5/3;
f.Position(4) = 1;

set(ax,"FontSize",8);

margin = struct("left",0.215,"right",.05,"top",.04,"bottom",.29);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

saveFigures(f,save_fig_opts)

%% Residuals of best pars by LL
[~,order] = sort(LL,"descend");

R = {zeros(nt,N,C.nsamps_per_condition),zeros(nt,numel(LL)- N,C.nsamps_per_condition)};
tt_min = E.t*1440;
C.ids = reshape(C.ids,[],C.nsamps_per_condition);
for i = 1:length(order)
    for j = 1:C.nsamps_per_condition
        load(sprintf("../../data/sims/%s/output_final.mat",C.ids(order(i),j)),"tracked")
        if i<=N
            R{1}(:,i,j) = interp1(round(tracked.t*1440),tracked.NT,tt_min)-E.D.A;
        else
            R{2}(:,i-N,j) = interp1(round(tracked.t*1440),tracked.NT,tt_min)-E.D.A;
        end
    end
end
R{1} = reshape(R{1},nt,[]);
R{2} = reshape(R{2},nt,[]);

%% make the figure of these residuals
f = figureOnRight("Name","TimeSeriesZScores_ByLikelihood_New");
ax = gobjects(1,nt-1);
for j = 1:nt-1
    tind = j+1;
    min_temp = min(min(R{1}(tind,:)),min(R{2}(tind,:))); % minimum residual at time tt(tind)
    max_temp = max(max(R{1}(tind,:)),max(R{2}(tind,:))); % maximum residual at time tt(tind)

    ax(j) = subplot(1,nt-1,j); hold on;
    [~,binEdges] = histcounts([R{1}(tind,:),R{2}(tind,:)]./E.D.S(tind));
    histogram(ax(j),R{1}(tind,:)./E.D.S(tind),"BinEdges",binEdges,"FaceColor",most_likely_color,"Normalization","pdf","DisplayName","Most Likely","EdgeColor","none")
    histogram(ax(j),R{2}(tind,:)./E.D.S(tind),"BinEdges",binEdges,"FaceColor",least_likely_color,"Normalization","pdf","DisplayName","Least Likely","EdgeColor","none")
    if min_temp~=0 || max_temp~=0
        xlim(max(abs([min_temp,max_temp])) * [-1 1]/E.D.S(tind)) % set xlim to be big enough to get largest |residual| and center at x=0
    end
    title(sprintf("t = %dh",round(24*E.t(tind))))
    xlabel("Z Score")
end
% legend(ax(2))
f.Units = "inches";
f.Position(3) = 5;
f.Position(4) = 1;
set(ax,"FontSize",8);

%%
margin = struct("left",0.02,"right",0.01,"top",.14,"bottom",0.29);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)

%% Residuals of accepted pars by LL
R = {zeros(nt,N,C.nsamps_per_condition),zeros(nt,numel(LL)- N,C.nsamps_per_condition)};
tt_min = E.t*1440;
C.ids = reshape(C.ids,[],C.nsamps_per_condition);
i_accepted = 1;
i_rejected = 1;
for i = 1:length(order)
    for j = 1:C.nsamps_per_condition
        load(sprintf("../../data/sims/%s/output_final.mat",C.ids(i,j)),"tracked")
        if accepted_parameters(i)
            R{1}(:,i_accepted,j) = interp1(round(tracked.t*1440),tracked.NT,tt_min)-E.D.A;
        else
            R{2}(:,i_rejected,j) = interp1(round(tracked.t*1440),tracked.NT,tt_min)-E.D.A;
        end
    end
    if accepted_parameters(i)
        i_accepted = i_accepted+1;
    else
        i_rejected = i_rejected+1;
    end
end
R{1} = reshape(R{1},nt,[]);
R{2} = reshape(R{2},nt,[]);

%% make the figure of these residuals
f = figureOnRight("Name","TimeSeriesZScores_ByAcceptance_New");
ax = gobjects(1,nt-1);
for j = 1:nt-1
    tind = j+1;
    min_temp = min(min(R{1}(tind,:)),min(R{2}(tind,:))); % minimum residual at time tt(tind)
    max_temp = max(max(R{1}(tind,:)),max(R{2}(tind,:))); % maximum residual at time tt(tind)

    ax(j) = subplot(1,nt-1,j); hold on;
    [~,binEdges] = histcounts([R{1}(tind,:),R{2}(tind,:)]./E.D.S(tind));
    histogram(ax(j),R{1}(tind,:)./E.D.S(tind),"BinEdges",binEdges,"FaceColor",accept_color,"Normalization","pdf","DisplayName","Accepted","EdgeColor","none")
    histogram(ax(j),R{2}(tind,:)./E.D.S(tind),"BinEdges",binEdges,"FaceColor",reject_color,"Normalization","pdf","DisplayName","Rejected","EdgeColor","none")
    if min_temp~=0 || max_temp~=0
        xlim(max(abs([min_temp,max_temp])) * [-1 1]/E.D.S(tind)) % set xlim to be big enough to get largest |residual| and center at x=0
    end
    title(sprintf("t = %dh",round(24*E.t(tind))))
    xlabel("Z Score","FontSize",8)
end
f.Units = "inches";
f.Position(3) = 4.5;
f.Position(4) = 1;

%%
margin = struct("left",0.03,"right",0.01,"top",.15,"bottom",0.3);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);
set(ax,"FontSize",8)

%%
saveFigures(f,save_fig_opts)

%% get abm trajectories
A = load(abm_data_file,"D");
accepted = A.D(:,accepted_parameters);
accepted = arrayify(accepted,"A",1);
accepted = sum(accepted,2);
rejected = A.D(:,~accepted_parameters);
rejected = arrayify(rejected,"A",1);
rejected = sum(rejected,2);
N = sum(accepted_parameters,"all");
[LL_sort,order] = sort(LL,"descend");
most_likely_inds = order(1:N);
most_likely = A.D(:,most_likely_inds);
most_likely = arrayify(most_likely,"A",1);
most_likely = sum(most_likely,2);
least_likely_inds = order(N+1:end);
least_likely = A.D(:,least_likely_inds);
least_likely = arrayify(least_likely,"A",1);
least_likely = sum(least_likely,2);

%% accepted vs most likely sample trajectories
acceptance_method = "all_profiles_resampled";

f = figureOnRight("Name","SamplesByAcceptedMostLikely_" + "New" + "_" + acceptance_method);
ax = gca; hold on;
nsamps = 40;

k = min(nsamps,N);
l_acc = plot(ax,E.t,squeeze(accepted(:,randperm(N,k))),"Color",[accept_color,nthroot(1/k,4)],"LineWidth",0.5);
k = min(nsamps,N);
l_mlik = plot(ax,E.t,squeeze(most_likely(:,randperm(N,k))),"Color",[most_likely_color,nthroot(1/k,4)],"LineWidth",0.5);
l_dat = plot(ax,E.t,E.D.A,"Color","black","LineStyle","--","LineWidth",0.5,"Marker","*","MarkerSize",3);

% legend(ax,[l_acc(1),l_mlik(1),l_dat],["Accepted","Rejected","Data"]);
set(ax,"FontSize",20)
xlabel(ax,"Time (d)")
ylabel(ax,"Cell Count")
set(ax,"FontSize",8)
f.Units = "inches";
f.Position(3) = 2.5;
f.Position(4) = 1;
ylim([0 1000])

%%
margin = struct("left",0.17,"right",.02,"top",.05,"bottom",0.3);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)

%% rejected vs least likely sample trajectories
acceptance_method = "all_profiles_resampled";

f = figureOnRight("Name","SamplesByRejectedLeastLikely_" + "New" + "_" + acceptance_method);
ax = gca; hold on;
nsamps = 40;

k = min(nsamps,length(least_likely_inds));
l_rej = plot(ax,E.t,squeeze(rejected(:,randperm(N,k))),"Color",[reject_color,nthroot(1/k,2)],"LineWidth",0.5);
k = min(nsamps,length(least_likely_inds));
l_least = plot(ax,E.t,squeeze(least_likely(:,randperm(length(least_likely_inds),k))),"Color",[least_likely_color,nthroot(1/k,2)],"LineWidth",0.5);
l_dat = plot(ax,E.t,E.D.A,"Color","black","LineStyle","--","LineWidth",0.5,"Marker","*","MarkerSize",3);

% legend(ax,[l_least(1),l_most(1),l_dat],["Accepted","Rejected","Data"]);
set(ax,"FontSize",20)
xlabel(ax,"Time (d)")
ylabel(ax,"Cell Count")
set(ax,"FontSize",8)
f.Units = "inches";
f.Position(3) = 2.5;
f.Position(4) = 1;
ylim([0 1000])

%%
margin = struct("left",0.17,"right",.02,"top",.05,"bottom",0.3);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)

%% least likely vs most likely sample trajectories
acceptance_method = "all_profiles_resampled";

f = figureOnRight("Name","SamplesByMostLeastLikely_" + "New" + "_" + acceptance_method);
ax = gca; hold on;
nsamps = 20;

k = min(nsamps,length(least_likely_inds));
l_least = plot(ax,E.t,squeeze(least_likely(:,randperm(length(least_likely_inds),k))),"Color",[least_likely_color,nthroot(1/k,4)],"LineWidth",0.5);
k = min(nsamps,N);
l_most = plot(ax,E.t,squeeze(most_likely(:,randperm(N,k))),"Color",[most_likely_color,nthroot(1/k,4)],"LineWidth",0.5);
l_dat = plot(ax,E.t,E.D.A,"Color","black","LineStyle","--","LineWidth",0.5,"Marker","*","MarkerSize",3);

% legend(ax,[l_least(1),l_most(1),l_dat],["Accepted","Rejected","Data"]);
set(ax,"FontSize",20)
xlabel(ax,"Time (d)")
ylabel(ax,"Cell Count")
set(ax,"FontSize",8)
f.Units = "inches";
f.Position(3) = 2;
f.Position(4) = 1;
ylim([0 1000])

%%
margin = struct("left",0.22,"right",.02,"top",.05,"bottom",0.3);
spacing = struct("horizontal",0.05,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%%
saveFigures(f,save_fig_opts)

%% set x data into log-lin scale
% x data is in log scale below the divider and linear scale above the divider
% the lbase is the logarithm used on the log scale side
function x = logLinScaleX(x,divider,lbase)

assert(isequal(x,sort(x))) % make sure x sorted data
assert(all(x>0)) % make sure all x>0 so the log plot makes sense

lowerX = x < divider;
if any(lowerX)
    maxLowerX = max(x(lowerX));
    x(lowerX) = log(x(lowerX))/log(lbase);
    c = divider - log(divider/maxLowerX)/log(lbase);
    x(lowerX) = x(lowerX) + c - max(x(lowerX));
end

end