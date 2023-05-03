% Actually going to draw from ABM parameter space based on SM parameters
% and SMoRe ParS profile likelihoods.

clearvars;

profile_file = "../ProfileLikelihood/data/MultiDimProfileLikelihoods.mat";
load(profile_file,"MDProfile","cohort_size","par_vals","profile_size")
