% This script sets up and calls the profile likelihood method for the ODE
% SM on a grid of ABM parameters. compare_every allows for comparing with a
% subset of ABM time points. This could be improved slightly by
% interpolating the ABM output at those times rather than the rounding I do
% to get tt. The function profileLikelihood computes the profile likelihood
% of all ODE model parameters at a given ABM parameter vector.

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../ODEFitting/OxControl/")

% folder to store plots and text files
cohort_name = "cohort_230124175743017";

% load data and compute the standard error
load("../../ODEFitting/OxControl/data/OptimalParameters_noapop.mat","P")
C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size");
Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"ode_state*");
load(sprintf("../../data/%s/output.mat",cohort_name),"ids");
load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");
t_abm = tracked.t;
compare_every = 6 / 24;

%% sample from the ABM output
tt = 0:compare_every:round(t_abm(end));
[~,tind] = min(abs(t_abm - tt),[],1);

Sum.ode_state_average = reshape(Sum.ode_state_average,numel(t_abm),2,[]);
Sum.ode_state_std = reshape(Sum.ode_state_std,numel(t_abm),2,[]);

%%
threshold = chi2inv(0.95,3); % compute threshold value for the parameter confidence intervals

% specify parameter ranges for bounds on profiling
para_ranges = [0,100;     % lambda
               0,100;  % alpha
               0,1e4];      % K

% set bounds for optimizing when profiling the other parameters
lb = [0;0;0];
ub = [Inf;Inf;1e4];
opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

% specify parameter names
para_names = {'\lambda', '\alpha', 'K'};
para_names_file_save = {'lambda', 'alpha', 'K'};

npars = size(P,1); % number of parameters

%%
Avg = Sum.ode_state_average(tind,:,:);
Std = Sum.ode_state_std(tind,:,:);
FF(1:size(Sum.ode_state_average,3)) = parallel.FevalFuture;
for i = 1:size(Sum.ode_state_average,3)
% for i = 3
    data = Avg(:,:,i);
    data_std = Std(:,:,i);
%     out(:,i) = profileLikelihood(P(:,i),tt,data,data_std,para_ranges,lb,ub,opts,threshold);
%     F = @(p) sum(((computeTimeSeries(p,tt) - data)./data_std).^2,'all');
    FF(i) = parfeval(@() profileLikelihood(P(:,i),tt,data,data_std,para_ranges,lb,ub,opts,threshold),1);
end

%%
out = cell(npars,size(Sum.ode_state_average,3));
for i = 1:size(Sum.ode_state_average,3)

    [idx,temp] = fetchNext(FF);
    out(:,idx) = temp;

    if mod(i,100)==0
        fprintf("Finished %d.\n",i)
    end

end

% save("ProfileLikelihoods.mat","out")

rmpath("../../ODEFitting/")

