clearvars
%% Program to run

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("..")

%% 2) Parameters : Please fill in
alpha = 0.05; % significance value for CI to determine if enough samples have been computed
options.limit_factor = 0.5; % how to set the limit for separating low and high impact factors
options.initialization_factor = 0.8; % how to determine when it has been sufficiently initialized
ci_relative_spread = 0.1; % how much the confidence interval can spread around the mean of the stochastic simulation output
nsim_max = 210;

nsamps = 10;
par_names = ["a","b"];
[D,I] = makeCMParameterDistributionsDictionary(par_names,distribution="Normal");

% Number of factors of uncertainty of the function studied :
nfac=numel(par_names); 

studied_function = @(x) moatSample(x,par_names,D,I,nsamps,alpha,ci_relative_spread);
[mu_star,sigma,order] = morris_simple(studied_function,nfac,1500);

%% reset path
rmpath("..")
