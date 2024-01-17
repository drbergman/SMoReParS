% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

% This is a new version that is being fit into the new workflow

clearvars;

addpath("../../../ODEFittingFns/")

file_name = "SMFitToData_New";

make_save = false;

optim_opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
fn = @computeTimeSeries;
fn_opts.condition_on_previous = false;
experimental_data = "data/ExperimentalData_New.mat";

%% set up inputs
[p,lb,ub] = basePars();
npars = length(p);
load(experimental_data,"t","D");
sm = struct("fn",fn,"opts",fn_opts);
F = @(p) getRawError(sm,p,t,D,[]);

%% optimize
[P,fstar] = fmincon(F,p,[],[],[],[],lb,ub,[],optim_opts);

%% save output
if make_save
    save("data/" + file_name,"P","fstar","fn_opts","lb","ub","fn","fn_opts","optim_opts") %#ok<*UNRCH>
end

%% reset path
rmpath("../../../ODEFittingFns/")
