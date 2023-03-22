clearvars;

% fit ABM parameters directly to data

cohort_name = "cohort_230124175743017";
C = load(sprintf("../data/%s/output.mat",cohort_name),"ids","cohort_size");

tt = [0      10     24     36     48     72    ]';        % hours
tt = tt/24;                                             % days

tt = round(1440*tt);
% Control data
data = [0.899  1.340  1.633  2.408  3.557  5.583]';   % millions of cells
data_std = [0.099  0.193  0.207  0.298  0.168  0.364]';   % millions of cells

factor = 100/data(1);

data_scaled = data * factor;
data_std_scaled = data_std * factor;

%% check with abm pars result in mean time series that pass with 1 SD of each data point; also compute mean RSS
C.ids = reshape(C.ids,prod(C.cohort_size),[]);
is_within_sd = false(size(C.ids,1),1); % abm pars that have mean within 1 SD
mean_residual = zeros(size(C.ids,1),1); % mean RSS of abm pars
% RSS = zeros(size(C.ids)); % RSS of each sample
for pi = 1:size(C.ids,1)
    par_time_series = zeros(length(tt),size(C.ids,2));
    for si = 1:size(C.ids,2)
    
        load(sprintf("../data/sims/%s/output_final.mat",C.ids(pi,si)),"tracked")
        par_time_series(:,si) = interp1(round(1440*tracked.t),tracked.NT,tt);

    end
    is_within_sd(pi) = all(abs(data_scaled-mean(par_time_series,2)) < data_std_scaled);
    mean_residual(pi) = mean(sum((abs(data_scaled-par_time_series)./data_std_scaled).^2,1));
    if mod(pi,250)==0
        fprintf("Finished %3.2f.\n",pi / size(C.ids,1))
    end

end

%% histogram of RSS by is_within_sd

figure; hold on;
hist_norm = "PDF";
histogram(mean_residual(is_within_sd),"FaceColor","green","normalization",hist_norm)
histogram(mean_residual(~is_within_sd),"FaceColor","red","normalization",hist_norm)
xlabel("RSS")
ylabel(hist_norm)

%% select by RSS cutoff
cutoff = 799/(3^7); % proportion of ABM parameters to select (the method using alphla, lambda, and K in a comination selected 799 of the 3^7 pars)
[RSS_sort,order] = sort(mean_residual);
v_cutoff = false(size(is_within_sd));
v_cutoff(order(1:round(cutoff*numel(order)))) = true;

%% histogram of RSS by cutoff

figure; hold on;
hist_norm = "count";
histogram(mean_residual(v_cutoff),"FaceColor","green","normalization",hist_norm)
histogram(mean_residual(~v_cutoff),"FaceColor","red","normalization",hist_norm)
xlabel("RSS")
ylabel(hist_norm)