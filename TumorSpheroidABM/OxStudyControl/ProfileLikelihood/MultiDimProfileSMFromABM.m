% This script sets up and calls the multi-dimensional profile likelihood method for the ODE
% SM on a grid of ABM parameters. compare_every allows for comparing with a
% subset of ABM time points. This could be improved slightly by
% interpolating the ABM output at those times rather than the rounding I do
% to get tt. The function profileLikelihood computes the profile likelihood
% of all ODE model parameters at a given ABM parameter vector.

clearvars;
addpath("../../../SurrogateModelFns/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions

cohort_name = "cohort_230124175743017";

files.data_file = sprintf("../../data/%s/summary_short.mat",cohort_name);
load(files.data_file,"cohort_size");

opts.force_serial = false;

n_sm_pars = 3;

n_vals_per_par = 10;

%% setup profile params
% specify parameter ranges for bounds on profiling
par_vals = {linspace(0.5,2.2,n_vals_per_par);     % lambda
               linspace(3.3,7.5,n_vals_per_par);  % alpha
               logspace(log10(300),log10(20000),n_vals_per_par)}; % K
par_names = ["\lambda","\alpha","K"];
profile_size = reshape(cellfun(@numel,par_vals),1,[]);
%% objfn_constants
objfn_constants.fn = @computeTimeSeries;
objfn_constants.fn_opts = [];
objfn_constants.weights = 1;
MDProfile = performMultiDimProfile(files,objfn_constants,par_vals,opts);

% save("data/MultiDimProfileLikelihoods.mat","MDProfile","par_vals","par_names","cohort_size","profile_size","-v7.3")

rmpath("../ODEFitting/")
rmpath("../../../SurrogateModelFns/")
rmpath("../../../ProfileLikelihoodFns/")


