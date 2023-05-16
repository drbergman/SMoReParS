% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;
profile_opts.force_serial = true;

addpath("../../ODEFittingFns/")


fn = @computeTimeSeries;
fn_opts.model_type = "logistic";

optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);


%% setup SM parameter and bounds
p = basePars(fn_opts.model_type);
switch fn_opts.model_type
    case "logistic"
        lb = [0;0];
        ub = [100;1e6];
    case "von_bertalanffy"
        lb = [0;1;0];
        ub = [100;Inf;100];
end
npars = numel(p);


%% load ABM data
data_file = "../PostAnalysis/data/summary.mat";

P = optimizeSMParsFromABM(data_file,p,fn,fn_opts,lb,ub,optim_opts,1,profile_opts);

save(sprintf("data/OptimalParameters_%s.mat",fn_opts.model_type),"P")

rmpath("../../ODEFittingFns/")
