% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../ODEFittingFns/")

par_names = ["\alpha","\theta","\beta"];

nsamps = 36;
par_file = "data/OptimalParameters.mat";
data_file = "../PostAnalysis/summary.mat";

fn = @computeTimeSeries;

f = testSMFitToABM(par_file,data_file,nsamps,fn,[],par_names);


