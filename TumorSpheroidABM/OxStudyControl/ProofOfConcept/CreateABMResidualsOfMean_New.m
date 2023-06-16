% This script will create the mean residuals of the ABM data compared to
% the experimental data


clearvars;
cohort_name = "cohort_230124175743017";

C = load(sprintf("../../data/%s/summary_short.mat",cohort_name));
E = load("../ODEFitting/data/ExperimentalData_New.mat","t","D");
nt = length(E.t);

n_abm_vecs = prod(C.cohort_size);
residuals_of_mean = zeros([nt,C.cohort_size]); % residuals of mean time series for each abm parameter vector
for pi = 1:numel(C.D) % (P)arameter (I)ndex
    residuals_of_mean(:,pi) = sum(C.D(pi).A,2) - E.D.A;
    if mod(pi,250)==0
        fprintf("Finished %3.2f.\n",pi / numel(C.D))
    end
end

%% save the output
save("data/ResidualsOfMean_New.mat","residuals_of_mean","-v7.3")
