% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;
force_serial = false;

addpath("../../ODEFittingFns/")

p = basePars();
npars = numel(p);

fn = @computeTimeSeries;

lb = zeros(3,1);
ub = [Inf;Inf;Inf];

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%% load ABM data
data_file = "../PostAnalysis/summary.mat";

P = optimizeSMParsFromABM(data_file,p,fn,[],lb,ub,opts,1,force_serial);

save("data/OptimalParameters.mat","P")

rmpath("../../ODEFittingFns/")
