% This script will use the combination of lambda and alpha found in the
% LambdaAlphaFn found from FitLambdaAlphaCurve.m and restricted by the 95%
% CI to restrict ABM parameter space. It will then sample points in ABM
% parameter space, classifying them as inside or outside the restricted
% region and plot their trajectories (if inside) and their RSS (both).

clearvars;

cohort_name = "cohort_230124175743017";
C = load(sprintf("../data/%s/output.mat",cohort_name),"cohort_size","lattice_parameters");
load("ProfileLikelihoods.mat")

npoints = size(out,2);
npars_ode = size(out,1);

%% create bounding hypersurfaces
BS = zeros(npars_ode,npoints,2);
for i = 1:npoints

    for j = 1:npars_ode
        BS(j,i,:) = out{j,i}(1,[1,end]);
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
load("ProfileLikelihoods_DataRestricted.mat","out");
load("LambdaAlphaFn.mat","f")
ODE_fit = load("../ODEFitting/ODEFittoData.mat");

% specify parameter ranges
para_ranges = [0,1e1;     % lambda
               0,1e1;  % alpha
               0,4e3];      % K

lb = [0;0;0];
ub = [1e1;1e1;4e3];
opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

threshold = chi2inv(0.95,3); % compute threshold value for the parameter confidence intervals

abm_region_1_log = false(size(Vq{1})); % bounded by lambda and alpha
for i = 1:size(out{1},2)
    temp = true(size(Vq{1}));
    for pi = 1:3
        temp = temp & Vq{pi,1} <= out{1}(pi,i) & Vq{pi,2} >= out{1}(pi,i);
    end
%     K_temp = profileLikelihood_one_ind(3,out{1}(1:3,i),ODE_fit.tt,ODE_fit.data,ODE_fit.data_std,para_ranges,lb,ub,opts,threshold,false);
%     plot(K_temp(1,:),K_temp(2,:))
%     drawnow
%     [~,min_ind] = min(K_temp(2,:),[],2);
%     temp = temp & Vq{3,1} <= K_temp(1,min_ind) & Vq{3,2} >= K_temp(1,min_ind);
    abm_region_1_log = abm_region_1_log | temp;
    disp(sum(abm_region_1_log,'all'))
end

% lambda = linspace(1,10,40000);
% lambda = [1.1,1.2,1.4,1.5341,1.6,2,2.5,3,5];
% load("LambdaAlphaFn.mat")
% alpha = f(lambda);
% lamalph = [lambda;alpha];
% for i = 1:numel(lambda)
%     temp = true(size(Vq{1}));
%     for pi = 1:2
%         temp = temp & Vq{pi,1} <= lamalph(pi,i) & Vq{pi,2} >= lamalph(pi,i);
%     end
%     abm_region_1_log = abm_region_1_log | temp;
% %     if mod(i,1000)==0
%         disp(sum(abm_region_1_log,'all'))
% %     end
% end

% load("ODEFittoData.mat","pstar"); % best ODE pars fit to data
% for pi = 1:2
%     abm_region_1_log = abm_region_1_log & Vq{pi,1} <= pstar(pi) & Vq{pi,2} >= pstar(pi);
% end
% 
% abm_region_2_log = abm_region_1_log & Vq{3,1} <= pstar(3) & Vq{3,2} >= pstar(3); % also bounded by hypersurface from K

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



