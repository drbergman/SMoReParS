% This script will compare fitting the ABM using SMoRe ParS against using
% the data directly.

clearvars;
addpath("~/Documents/MATLAB/myfunctions/")
cohort_name = "cohort_230124175743017";

%% load SMoRe ParS-derived fitting
load("../../ProfileLikelihood/OxControl/data/ABMParamEstimates_FromProfile_WithK.mat","LP1","abm_region_1_log")

%% load cohort data and experimental data
C = load(sprintf("../../data/%s/output.mat",cohort_name),"ids","cohort_size");
load("../../ODEFitting/OxControl/data/ExperimentalData.mat")
n = size(count,1);
tt = round(tt*1440); % time in minutes to avoid rounding errors
nt = length(tt);

%% load residuals of mean
load("data/ResidualsOfMean.mat")

%% Compute log-likelihood values
LL = -0.5*(n*log(2*pi()) + sum(count_std.^2) + sum((residuals_of_mean./count_std).^2,1));

%% stratify LL by SMoRe ParS Acceptance/Rejection
LL_accepted = LL(abm_region_1_log(:));
LL_rejected = LL(~abm_region_1_log(:));

%% plot histograms of both
figure; hold on;
ax = gca;
[~,binEdges] = histcounts(LL);
histogram(LL_accepted,"BinEdges",binEdges,"FaceColor","green","Normalization","pdf","DisplayName","Selected Samples","EdgeColor","none")
histogram(LL_rejected,"BinEdges",binEdges,"FaceColor","red","Normalization","pdf","DisplayName","Rejected Samples","EdgeColor","none")
legend(ax,"Location","northwest","FontSize",16)
xlabel("log-likelihood","FontSize",16)
ylabel("PDF","FontSize",16)

%% distribution relative likelihood to best setup
N = sum(abm_region_1_log,"all");
divider = 1e-2;
minRelProb = exp(min(LL)-max(LL));
lbase = divider/minRelProb;
fn = @logLinScaleX;
XData = {exp(sort(LL,"ascend")-max(LL)),exp(sort(LL_accepted)-max(LL)),exp(sort(LL_rejected)-max(LL))};
CDF = {linspace(0,1,numel(LL)),linspace(0,1,length(LL_accepted)),linspace(0,1,length(LL_rejected))};
Names = ["All","Accepted","Rejected"];
Colors = ["b","g","r"];

%% plot cdf
figure;
hold on;
for i = 1:3
plot(fn(XData{i},divider,lbase),CDF{i},"DisplayName",Names(i),"Color",Colors(i),"LineWidth",2)
end
xlabel("Relative Likelihood","FontSize",16)
ylabel("CDF","FontSize",16)
legend(gca,"location","best","FontSize",16,"AutoUpdate","off")
xL = xlim;
xL(1) = fn(minRelProb,divider,lbase);
minPow10 = log10(minRelProb);
maxPow10 = log10(divider);
lowerPow10 = flip(maxPow10:ceil(0.5*(minPow10-maxPow10)):minPow10);
lowerXTicks = fn(10.^lowerPow10,divider,lbase);
upperXTicks = .2:.2:1;
xticks(unique([lowerXTicks,upperXTicks]))
xticklabels([divider*lbase.^(lowerXTicks-divider),upperXTicks])
xL(2) = 1;
xlim(xL)
xline(divider,"LineStyle","--")
annotation("textarrow",[.2875,.1875],[0.9,0.9],"String","Log Scale");
annotation("arrow",[0.3675,0.4675],[0.9,0.9]);
annotation("textarrow",[.7625,.8625],[0.15,0.15],"String","Linear Scale");
annotation("arrow",[0.6675,0.5675],[0.15,0.15]);
savefig("figs/RelativeLikelihoodCDF")
print("figs/RelativeLikelihoodCDF","-dpng")

%% attempt to plot PDF of above
figure;
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
figure;
hold on;
% x = [0,logspace(-30,0,20)];
x = linspace(0,1,26);
for i = 1:3
    temp = histcounts(XData{i},x,"Normalization","pdf");
    plot(.5*(x(1:end-1)+x(2:end)),temp,"DisplayName",Names(i),"Color",Colors(i),"LineWidth",2)
end
set(gca,"YScale","log")
xlabel("Rleative Likelihood","FontSize",16)
ylabel("PDF","FontSize",16)
legend(gca,"location","best","FontSize",16,"AutoUpdate","off")
savefig("figs/RelativeLikelihoodPDF")
print("figs/RelativeLikelihoodPDF","-dpng")

%% plotting percentage of N best LLs that were accepted
figure;
[~,order] = sort(LL,"descend");
accepted_sort = abm_region_1_log(order);
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

savefig("figs/AcceptanceOfMostLikely")
print("figs/AcceptanceOfMostLikely","-dpng")

%% reorient above to have categories be most likely and least likely and compare 
figure;
ax=gca;
[LL_sort,order] = sort(LL,"ascend");
accepted_sort = abm_region_1_log(order);
n_quants = 47;
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
bar(linspace(0,100,n_quants),100*y,"stacked","BarWidth",1,"EdgeColor","none");
xlim([0,100]+50/(n_quants-1)*[-1 1])
ylim([0 100])
legend(["Accepted","Rejected"],"location","northwest","AutoUpdate","off")
ylabel("Probability")
xlabel("Likelihood Percentile")
set(ax,"FontSize",16)
yticks(ax,0:25:100)
ax.XAxis.TickLength = [0 0];
xline(100*(1-N/numel(LL)))

savefig("figs/AcceptanceByLikelihoodQuantile")
print("figs/AcceptanceByLikelihoodQuantile","-dpng")

%% Residuals of best pars by LL
[~,order] = sort(LL,"descend");

R = {zeros(length(tt),N),zeros(length(tt),numel(LL)-N)};
tt_min = tt;
for i = 1:length(order)
    load(sprintf("../../data/sims/%s/output_final.mat",C.ids(order(i))),"tracked")
    if i<=N
        R{1}(:,i) = interp1(round(tracked.t*1440),tracked.NT,tt_min)-count;
    else
        R{2}(:,i-N) = interp1(round(tracked.t*1440),tracked.NT,tt_min)-count;
    end
end

%% make the figure of these residuals
figure;
ax = gobjects(length(tt),1);
for j = 1:length(tt)
    min_temp = min(min(R{1}(j,:)),min(R{2}(j,:))); % minimum residual at time tt(j)
    max_temp = max(max(R{1}(j,:)),max(R{2}(j,:))); % maximum residual at time tt(j)

    ax(j) = subplot(1,length(tt),j); hold on;
    [~,binEdges] = histcounts([R{1}(j,:),R{2}(j,:)]./count_std(j));
    histogram(ax(j),R{1}(j,:)./count_std(j),"BinEdges",binEdges,"FaceColor","green","Normalization","pdf","DisplayName","Most Likely","EdgeColor","none")
    histogram(ax(j),R{2}(j,:)./count_std(j),"BinEdges",binEdges,"FaceColor","red","Normalization","pdf","DisplayName","Least Likely","EdgeColor","none")
    if min_temp~=0 || max_temp~=0
        xlim(max(abs([min_temp,max_temp])) * [-1 1]/count_std(j)) % set xlim to be big enough to get largest |residual| and center at x=0
    end
    title(sprintf("t = %3.2f d",tt(j)/1440))
    xlabel("Z Score")
end
legend(ax(2))

savefig("figs/TimeSeriesZScores_ByLikelihood")
print("figs/TimeSeriesZScores_ByLikelihood","-dpng")

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