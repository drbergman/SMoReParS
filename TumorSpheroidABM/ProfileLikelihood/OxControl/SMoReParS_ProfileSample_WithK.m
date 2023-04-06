% This script will use the combination of lambda and alpha found in the
% LambdaAlphaFn found from FitLambdaAlphaCurve.m and restricted by the 95%
% CI to restrict ABM parameter space. It will then sample points in ABM
% parameter space, classifying them as inside or outside the restricted
% region and plot their trajectories (if inside) and their RSS (both).

clearvars;

cohort_name = "cohort_230124175743017";
C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size","lattice_parameters");
load("data/ProfileLikelihoods.mat")

npoints = size(out,2);
npars_ode = size(out,1);

%% create bounding hypersurfaces
BS = zeros(npars_ode,npoints,2);
threshold = chi2inv(0.95,3);
for i = 1:npoints
    for j = 1:npars_ode
        [BS(j,i,1),BS(j,i,2)] = getProfileBounds(out{j,i},threshold);
    end
end

BS = reshape(BS,[npars_ode,C.cohort_size,2]);

%% interpolate at finer grid

pars = {C.lattice_parameters.values};
npars_abm = numel(pars);
nx = 3; % probably want to choose odd so the 3 calculated values along each dimension are used
par_grid = cell(npars_abm,1);
for i = 1:npars_abm
    par_grid{i} = linspace(pars{i}(1),pars{i}(end),nx);
end

PG = cell(npars_abm,1);
[PG{:}] = ndgrid(par_grid{:});

Vq = cell(npars_ode,2);
for pi = 1:npars_ode
    Vq{pi,1} = interpn(pars{:},squeeze(BS(pi,:,:,:,:,:,:,:,1)),PG{:});
    Vq{pi,2} = interpn(pars{:},squeeze(BS(pi,:,:,:,:,:,:,:,2)),PG{:});
end

%% see which abm pars have the best ode par within the bounding hypersurfaces
PLFromData = load("data/ProfileLikelihoods_DataRestricted.mat","out");
% load("data/LambdaAlphaFn.mat","f")
% ODE_fit = load("../../ODEFitting/OxControl/data/ODEFittoData.mat");

% specify parameter ranges
para_ranges = [0,1e1;     % lambda
               0,1e1;  % alpha
               0,4e3];      % K

lb = [0;0;0];
ub = [1e1;1e1;4e3];
opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

abm_region_1_log = false(size(Vq{1})); % bounded by lambda, alpha, and K
[lambda_start,lambda_end] = getProfileBounds(PLFromData.out{1}([1,end],:),threshold);
% lambda_start = PLFromData.out{1}(1,1);
% lambda_end = PLFromData.out{1}(1,end);
for i = 1:size(PLFromData.out{1},2)
    if PLFromData.out{1}(1,i)<lambda_start || PLFromData.out{1}(1,i)>lambda_end
        continue;
    end
    temp = true(size(Vq{1}));
    for pi = 1:3
        temp = temp & Vq{pi,1} <= PLFromData.out{1}(pi,i) & Vq{pi,2} >= PLFromData.out{1}(pi,i);
    end
    abm_region_1_log = abm_region_1_log | temp;
end

%% store the parameter vectors by both restrictions
I = cell(7,1);
[I{:}] = ind2sub(size(abm_region_1_log),find(abm_region_1_log));

LP1 =  C.lattice_parameters;
for i = 1:7
    LP1(i).values = reshape(par_grid{i}(I{i}),[],1);
end

%%
if nx~=3
    warning("abm_region_1_log does not correspond to the sampled 3^7 grid.")
end
save("ABMParamEstimates_FromProfile_WithK","LP1","abm_region_1_log")



