% Actually going to draw from ABM parameter space based on SM parameters
% and SMoRe ParS profile likelihoods.

clearvars;
addpath("../../../ProfileLikelihoodFns/")

%% load ABM profile data
profile_file = "../ProfileLikelihood/data/MultiDimProfileLikelihoods.mat";
load(profile_file,"MDProfile")

%% load rest of profile data
load(profile_file,"cohort_size","par_vals","profile_size")

%% load cohort data
cohort_name = "cohort_230124175743017";
CC=load(sprintf("../../data/%s/summary_short.mat",cohort_name),"vals");
abm_par_vals = CC.vals;
npars_abm = length(abm_par_vals);
clear CC

%% load SM profile data
PLFromData = load("../ProfileLikelihood/data/ProfileLikelihoods_DataRestricted.mat","out");

%% sample SM pars using likelihood
ll_profile_irregular_partition = -0.5*PLFromData.out{1}(end,:);
lambda_grid = linspace(PLFromData.out{1}(1,1),PLFromData.out{1}(1,end),1001); % obviate need for integrating by getting a uniform grid
ll_profile = interp1(PLFromData.out{1}(1,:),ll_profile_irregular_partition,lambda_grid); % get ll values on grid
profile_cdf = cumsum(exp(ll_profile-max(ll_profile))); % convert to likelihood values (relative to best) and "integrate"
profile_cdf = profile_cdf/profile_cdf(end); % normalize to a true CDF
n_sm_samps = 1000;
u = rand(n_sm_samps,1);
lambda_sample = interp1(profile_cdf,lambda_grid,u);
alpha_sample = interp1(PLFromData.out{1}(1,:),PLFromData.out{1}(2,:),lambda_sample);
K_sample = interp1(PLFromData.out{1}(1,:),PLFromData.out{1}(3,:),lambda_sample);

%% post-data loading stuff
npars_sm = length(profile_size);
n_abm_vecs = prod(cohort_size);
MDProfile = reshape(MDProfile,[prod(cohort_size),profile_size]);
true_par_vals = [lambda_sample,alpha_sample,K_sample];
bad_samps = lambda_sample < par_vals{1}(1) | lambda_sample > par_vals{1}(end) ...
    | alpha_sample < par_vals{2}(1) | alpha_sample > par_vals{2}(end) ...
    | K_sample < par_vals{3}(1) | K_sample > par_vals{3}(end);
true_par_vals(bad_samps,:) = [];
n_sm_samps = size(true_par_vals,1);

%% interpolate
LL_est = zeros(n_abm_vecs,n_sm_samps);
for i = 1:n_sm_samps
    LL_est(:,i) = interpn(1:2187,par_vals{:},MDProfile,1:2187,true_par_vals(i,1),true_par_vals(i,2),true_par_vals(i,3));
end
figure;
S = sort(LL_est,1);
P = exp(S-S(end,:));
plot(S-S(end,:),linspace(0,1,size(S,1))',"LineWidth",0.25,"Marker",".")

%% pick ABM grid point to be base point
abm_grid_ind = zeros(n_sm_samps,1);
for i = 1:n_sm_samps
    abm_grid_ind(i) = randsample(n_abm_vecs,1,true,exp(LL_est(:,i)-S(end,i)));
end

%% pick direction to search in for each parameter
LL_est = reshape(LL_est,cohort_size);
X=cell(npars_abm,1);
Y = cell(npars_abm,1);
[X{:}] = ind2sub(cohort_size,abm_grid_ind);
% par_dir = zeros(npars_abm,1);
for i = 1:npars_abm
    Y{i} = X{i}*ones(npars_abm,1);
    if X{i}==1
        Y{i}(i) = 2;
    elseif X{i} == cohort_size(i)
        Y{i}(i) = X{i}-1;
    else
        temp = LL_est;
        for j = 1:npars_abm
            if j~=i
                temp = sliceof(temp,j,X{j});
            end
        end
        ll_neighbors = temp(X{i}+[-1,1]);
        ll_neighbors = ll_neighbors - max(ll_neighbors);
        Y{i}(i) = randsample(X{i}+[-1,1],1,true,exp(ll_neighbors));
    end
end

%% draw each parameter from this hyperrectangle
neighbor_inds = sub2ind(cohort_size,Y{:});
G = LL_est(neighbor_inds) - LL_est(abm_grid_ind);
ns = 5;
c = log(1+(exp(G)-1).*rand(npars_abm,ns))./G;
sampled_abm_par_vals = zeros(ns,npars_abm);
for i = 1:npars_abm
    sampled_abm_par_vals(:,i) = interp1([0 1],abm_par_vals{i}([X{i},Y{i}(i)]),c(i,:));
end

%% plot cdf of each par
figureOnRight;
for i = 1:npars_abm
    subplot(npars_abm,1,i);
    plot(sort(sampled_abm_par_vals(:,i)),linspace(0,1,ns));
end


