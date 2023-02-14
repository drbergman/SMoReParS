clearvars;

% fit ABM parameters directly to data

cohort_name = "cohort_230124175743017";
C = load(sprintf("../data/%s/output.mat",cohort_name),"ids");

tt = [0      10     24     36     48     72    ]';        % hours
tt = tt/24;                                             % days

tt = round(1440*tt);
% Control data
data = [0.899  1.340  1.633  2.408  3.557  5.583]';   % millions of cells
data_std = [0.099  0.193  0.207  0.298  0.168  0.364]';   % millions of cells

factor = 100/data(1);

data_scaled = data * factor;
data_std_scaled = data_std * factor;

v = zeros(size(C.ids));
for i = 1:numel(C.ids)
    
    load(sprintf("../data/sims/%s/output_final.mat",C.ids(i)),"tracked")
    temp = interp1(round(1440*tracked.t),tracked.NT,tt);
    v(i) = all(abs(data_scaled-temp) < data_std_scaled);

    if mod(i,1000)==0
        fprintf("Finished %3.2f.\n",i / numel(C.ids))
    end

end

%% residuals

res = zeros(length(tt),numel(C.ids));
for i = 1:numel(C.ids)
    
    load(sprintf("../data/sims/%s/output_final.mat",C.ids(i)),"tracked")
    temp = interp1(round(1440*tracked.t),tracked.NT,tt);
    res(:,i) = abs(data_scaled-temp) ./ data_std_scaled;

    if mod(i,1000)==0
        fprintf("Finished %3.2f.\n",i / numel(C.ids))
    end

end

%% RSS

RSS = sum(res.^2,1);

%% histogram of RSS by v

v = v==1;
figure; hold on;
hist_norm = "PDF";
histogram(RSS(v),"FaceColor","green","normalization",hist_norm)
histogram(RSS(~v),"FaceColor","red","normalization",hist_norm)
xlabel("RSS")
ylabel(hist_norm)

%% select by RSS cutoff
cutoff = 0.1; % proportion of ABM parameters to select
[RSS_sort,order] = sort(RSS);
v_cutoff = false(size(v));
v_cutoff(order(1:round(cutoff*numel(order)))) = true;

%% histogram of RSS by cutoff

figure; hold on;
hist_norm = "PDF";
histogram(RSS(v_cutoff),"FaceColor","green","normalization",hist_norm)
histogram(RSS(~v_cutoff),"FaceColor","red","normalization",hist_norm)
xlabel("RSS")
ylabel(hist_norm)