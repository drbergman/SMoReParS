% This is a first attempt to draw values from ABM parameter space using the
% full likelihood (what I'm calling here the multi-dimensional profile or
% MDProfile)

clearvars;

profile_file = "../ProfileLikelihood/data/MultiDimProfileLikelihoods.mat";

load(profile_file,"MDProfile","cohort_size","par_vals","profile_size")

%% load cohort data
cohort_name = "cohort_230124175743017";
CC=load(sprintf("../../data/%s/summary_short.mat",cohort_name),"vals");
abm_par_vals = CC.vals;
npars_abm = length(abm_par_vals);
clear CC
%% post-data loading stuff

npars_sm = length(profile_size);
n_abm_vecs = prod(cohort_size);
MDProfile = reshape(MDProfile,[prod(cohort_size),profile_size]);
true_par_vals = zeros(npars_sm,1);
for i = 1:npars_sm
    interval_length = par_vals{i}(end)-par_vals{i}(1);
    true_par_vals(i) = interval_length*rand() + par_vals{i}(1);
end

%% interpolate
LL_est = interpn(1:2187,par_vals{:},MDProfile,1:2187,true_par_vals(1),true_par_vals(2),true_par_vals(3));
figure;
S = sort(LL_est);
P = exp(S-S(end));
plot(P,linspace(0,1,numel(S)),"LineWidth",0.25,"Marker",".")

%% find reasonable neighborhood
LL_est = reshape(LL_est,cohort_size);
LL_rel = LL_est - max(LL_est,[],"all");
NZ = find(exp(LL_rel)>0); % grid points with nonzero probability
X = cell(npars_abm,1);
refine_factor = 5;
for i = 1:numel(NZ)
    LL_rel_ref = LL_rel;
    abm_par_vals_refined = abm_par_vals;
    while true
        [~,temp_ind] = max(LL_rel_ref,[],"all");
        [X{:}] = ind2sub(size(LL_rel_ref),temp_ind);
        abm_par_grid_inds = cell(npars_abm,1);
        abm_par_grid_vals = cell(npars_abm,1);
        for j = 1:npars_abm
            abm_par_grid_inds{j} = unique(min(length(abm_par_vals_refined{j}),max(1,X{j}+(-1:1))));
            abm_par_grid_vals{j} = abm_par_vals_refined{j}(abm_par_grid_inds{j});
            abm_par_vals_refined{j} = linspace(abm_par_grid_vals{j}(1),abm_par_grid_vals{j}(end),refine_factor);
        end
        [X{:}] = ndgrid(abm_par_vals_refined{:});
        LL_rel_ref = interpn(abm_par_grid_vals{:},LL_rel_ref(abm_par_grid_inds{:}),X{:});
        if sum(exp(LL_rel_ref)>0,"all")>1
            break;
        end
    end
end


%% 
