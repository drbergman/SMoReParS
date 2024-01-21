% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;
force_serial = true;
raw_error_opts.resample = true;
raw_error_opts.t = 15:15:75;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../SurrogateModelFns/")

sm.opts.model_type = "exponential";
sm.custom_solve_sm_fn = @customSolveSM;

optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);


%% setup SM parameter and bounds
p = basePars(sm.opts.model_type);
switch sm.opts.model_type
    case "exponential"
        lb = 0;
        ub = 0.2; % y(75) = y0*exp(lambda * 75) = y0*(3.7e32)^lambda grows pretty fast in lambda!
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
    resample_t=15:15:75);

save(sprintf("data/OptimalParameters_%s.mat",sm.opts.model_type),"P","sm")

rmpath("../../SurrogateModelFns/")
