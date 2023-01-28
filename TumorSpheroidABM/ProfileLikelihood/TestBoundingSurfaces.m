% I used this script to help find a way to sample from ABM parameter space
% given the hypersurfaces. First, create the bounding hypersurfaces. Then,
% interpolate them at a fine grid (optional). 

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
load("../ODEFitting/ODEFittoData.mat","pstar")

Z1 = true(size(Vq{1})); % bounded by lambda and alpha
for pi = 1:2
    Z1 = Z1 & Vq{pi,1} <= pstar(pi) & Vq{pi,2} >= pstar(pi);
end

Z2 = Z1 & Vq{3,1} <= pstar(3) & Vq{3,2} >= pstar(3); % also bounded by hypersurface from K

%% store the parameter vectors by both restrictions
I = cell(7,1);
[I{:}] = ind2sub(size(Z1),find(Z1));

I2 = cell(7,1);
[I2{:}] = ind2sub(size(Z2),find(Z2));

LP1 =  C.lattice_parameters;
for i = 1:7
    LP1(i).values = reshape(par_grid{i}(I{i}),[],1);
end

LP2 =  C.lattice_parameters;
for i = 1:7
    LP2(i).values = reshape(par_grid{i}(I2{i}),[],1);
end

%%
save("ABMParamEstimates","LP1","LP2")


