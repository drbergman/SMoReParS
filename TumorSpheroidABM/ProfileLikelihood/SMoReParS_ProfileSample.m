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
nx = 5; % probably want to choose odd so the 3 calculated values along each dimension are used
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

Z1 = false(size(Vq{1})); % bounded by lambda and alpha
for i = 1:size(out{1},2)
    temp = true(size(Vq{1}));
    for pi = 1:2
        temp = temp & Vq{pi,1} <= out{1}(pi,i) & Vq{pi,2} >= out{1}(pi,i);
    end
    Z1 = Z1 | temp;
%     disp(sum(Z1,'all'))
end
Z1 = Z1 & Vq{3,1} >= 597; % the K profile says don't go below 597

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
%     Z1 = Z1 | temp;
% %     if mod(i,1000)==0
%         disp(sum(Z1,'all'))
% %     end
% end
% for pi = 1:2
%     Z1 = Z1 & Vq{pi,1} <= pstar(pi) & Vq{pi,2} >= pstar(pi);
% end
% 
% Z2 = Z1 & Vq{3,1} <= pstar(3) & Vq{3,2} >= pstar(3); % also bounded by hypersurface from K

%% store the parameter vectors by both restrictions
I = cell(7,1);
[I{:}] = ind2sub(size(Z1),find(Z1));

LP1 =  C.lattice_parameters;
for i = 1:7
    LP1(i).values = reshape(par_grid{i}(I{i}),[],1);
end

%%
save("ABMParamEstimates_FromProfile2","LP1","Z1")


