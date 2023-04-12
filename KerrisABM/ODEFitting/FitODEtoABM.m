% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;
force_serial = false;

addpath("../../ODEFittingFns/")

p = basePars();
npars = numel(p);

fn = @computeTimeSeries;

lb = zeros(3,1);
ub = [Inf;1;Inf]; % even if I let theta go to Inf, values bigger than 1 are not selected (perhaps they could be if alpha << 1? but that's not really a parameter range for theta we are interested in)

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%% load ABM data
data_file = "../PostAnalysis/summary.mat";

P = optimizeSMParsFromABM(data_file,p,fn,[],lb,ub,opts,1,force_serial);

save("data/OptimalParameters.mat","P")

rmpath("../../ODEFittingFns/")
