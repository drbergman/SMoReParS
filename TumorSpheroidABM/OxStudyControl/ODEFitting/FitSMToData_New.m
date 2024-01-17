% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

% This is a new version that is being fit into the new workflow
clear all
clearvars;
addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../SurrogateModelFns/")

%% reset persistent variables first
clear rawError solveSM

%% continue...
file_name = "SMFitToData_New";
experimental_data = "data/ExperimentalData_New.mat";

%% set up inputs
[p,lb,ub] = basePars();
npars = length(p);
load(experimental_data,"t","D");

make_save = false;

optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

% Using base method
% sm.type = "ode";
% sm.solver = @ode45;
% sm.fn = @odefn;
% sm.t0 = 0;
% sm.y0 = [90;10];
% sm.post_processor = @(x) sum(deval(t,x),1)';

% Using custom solve sm method
% sm.custom_solve_sm = @customSolveSM;

% Using custom raw error method
sm.custom_raw_error_fn = @customRawError;

if isfield(sm,"custom_raw_error_fn")
    F = @(p) sm.custom_raw_error_fn(sm,p,t,D,[],struct());
else
    F = @(p) getRawError(sm,p,t,D,[]);
end

%% optimize
[P,fstar] = fmincon(F,p,[],[],[],[],lb,ub,[],optim_opts);

%% save output
if make_save
    save("data/" + file_name,"P","fstar","fn_opts","lb","ub","fn","fn_opts","optim_opts") %#ok<*UNRCH>
end

%% reset path
rmpath("../../../SurrogateModelFns/")
