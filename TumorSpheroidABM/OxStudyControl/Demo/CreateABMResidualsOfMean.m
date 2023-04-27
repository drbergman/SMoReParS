% This script will create the mean residuals of the ABM data compared to
% the experimental data


clearvars;
cohort_name = "cohort_230124175743017";

C = load(sprintf("../../data/%s/output.mat",cohort_name),"ids","cohort_size");
load("../../OxStudyControl/ODEFitting/data/ExperimentalData.mat","tt","count")
tt = round(tt*1440); % time in minutes to avoid rounding errors
nt = length(tt);

C.ids = reshape(C.ids,prod(C.cohort_size),[]);
residuals_of_mean = zeros(nt,size(C.ids,1)); % residuals of mean time series for each abm parameter vector
for pi = 1:size(C.ids,1) % (P)arameter (I)ndex
    par_time_series = zeros(length(tt),size(C.ids,2));
    for si = 1:size(C.ids,2)
        load(sprintf("../../data/sims/%s/output_final.mat",C.ids(pi,si)),"tracked")
        par_time_series(:,si) = interp1(round(1440*tracked.t),tracked.NT,tt);
    end
    par_time_series = mean(par_time_series,2);
    residuals_of_mean(:,pi) = par_time_series-count;
    if mod(pi,250)==0
        fprintf("Finished %3.2f.\n",pi / size(C.ids,1))
    end
end

%% save the output
save("data/ResidualsOfMean.mat","residuals_of_mean","-v7.3")
