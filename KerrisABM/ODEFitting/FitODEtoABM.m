% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;
force_serial = true;
raw_error_opts.resample = true;
raw_error_opts.t = 15:15:75;

addpath("../../ODEFittingFns/")


sm.fn = @computeTimeSeries;
sm.opts.model_type = "logistic";

optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);


%% setup SM parameter and bounds
p = basePars(fn_opts.model_type);
switch fn_opts.model_type
    case "logistic"
        lb = [0;0];
        ub = [1;1e6];
    case "von_bertalanffy"
        lb = [0;1;0];
        ub = [100;Inf;100];
end
npars = numel(p);


%% load ABM data
files.data = "../PostAnalysis/data/summary.mat";

P = optimizeSMParsFromABM(files,sm,p,lb,ub,optim_opts,1,force_serial=force_serial,...
    raw_error_opts=raw_error_opts);

save(sprintf("data/OptimalParameters_%s.mat",fn_opts.model_type),"P")

rmpath("../../ODEFittingFns/")
