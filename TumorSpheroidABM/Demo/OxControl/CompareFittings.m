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

