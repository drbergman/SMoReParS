% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;
force_serial = true;
cohort_name = "cohort_230124175743017";

addpath("../../../ODEFittingFns/")

p = basePars();
npars = numel(p);

fn = @computeTimeSeries;

lb = [0;0;0];
ub = [Inf;Inf;1e4];

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%% load ABM data
data_file = sprintf("../../data/%s/summary.mat",cohort_name);

P = optimizeSMParsFromABM(data_file,p,fn,[],lb,ub,opts,1,force_serial);

% save("data/OptimalParameters_Using_optimizeSMParsFromABM.mat","P")

rmpath("../../../ODEFittingFns/")
