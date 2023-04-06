clearvars;

% fit ABM parameters directly to data

cohort_name = "cohort_230124175743017";
C = load(sprintf("../../data/%s/output.mat",cohort_name),"ids","cohort_size");

load("data/ExperimentalData.mat")
tt = round(tt*1440); % time in minutes to avoid rounding errors
%% check with abm pars result in mean time series that pass with 1 SD of each data point; also compute mean log-likelihood (LL) 
C.ids = reshape(C.ids,prod(C.cohort_size),[]);
is_within_sd = false(size(C.ids,1),1); % abm pars that have mean within 1 SD
mean_residual = zeros(size(C.ids,1),1); % mean LL of abm pars
par_residual = zeros(size(C.ids,1),1); % mean LL of abm pars
% LL = zeros(size(C.ids)); % LL of each sample
for pi = 1:size(C.ids,1) % (P)arameter (I)ndex
    par_time_series = zeros(length(tt),size(C.ids,2));
    for si = 1:size(C.ids,2)
    
        load(sprintf("../../data/sims/%s/output_final.mat",C.ids(pi,si)),"tracked")
        par_time_series(:,si) = interp1(round(1440*tracked.t),tracked.NT,tt);

    end
    is_within_sd(pi) = all(abs(count-mean(par_time_series,2)) < count_std);
    mean_residual(pi) = mean(sum((abs(count-par_time_series)./count_std).^2,1));
    par_residual(pi) = sum(((count-mean(par_time_series,2))./count_std).^2,1);
    if mod(pi,250)==0
        fprintf("Finished %3.2f.\n",pi / size(C.ids,1))
    end

end

%% histogram of log-likelihood by is_within_sd

figure; hold on;
hist_norm = "PDF";
histogram(mean_residual(is_within_sd),"FaceColor","green","normalization",hist_norm)
histogram(mean_residual(~is_within_sd),"FaceColor","red","normalization",hist_norm)
xlabel("log-likelihood")
ylabel(hist_norm)

%% select by LL cutoff comparing to data first then taking mean across sims
cutoff = 799/(3^7); % proportion of ABM parameters to select (the method using alphla, lambda, and K in a comination selected 799 of the 3^7 pars)
[~,order] = sort(mean_residual);
v_cutoff = false(size(is_within_sd));
v_cutoff(order(1:round(cutoff*numel(order)))) = true;

%% histogram of LL by cutoff

figure; hold on;
hist_norm = "count";
histogram(mean_residual(v_cutoff),"FaceColor","green","normalization",hist_norm)
histogram(mean_residual(~v_cutoff),"FaceColor","red","normalization",hist_norm)
xlabel("log-likelihood")
ylabel(hist_norm)

%% select by LL cutoff taking mean of ABM sims first then comparing to data
cutoff = 799/(3^7); % proportion of ABM parameters to select (the method using alphla, lambda, and K in a comination selected 799 of the 3^7 pars)
[~,order] = sort(par_residual);
v_cutoff = false(size(is_within_sd));
v_cutoff(order(1:round(cutoff*numel(order)))) = true;

%% histogram of LL by cutoff
ll_const = numel(count)*log(2*pi()) + sum(log(count_std.^2));
data = -0.5*(ll_const+par_residual);
figure; hold on;
hist_norm = "count";
histogram(data(v_cutoff),"FaceColor","green","normalization",hist_norm)
histogram(data(~v_cutoff),"FaceColor","red","normalization",hist_norm)
xlabel("log-likelihood")
ylabel(hist_norm)