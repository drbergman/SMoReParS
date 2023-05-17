% this script identify admissible points in ABM parameter space from SMoRe 
% ParS.

clearvars;

addpath("../../../ProfileLikelihoodFns/")
addpath("../../../SampleABMFns/")

cohort_name = "cohort_2303301105";

files.abm_data_file = sprintf("../../data/%s/summary_new.mat",cohort_name);
files.profile_from_abm_file = "../ProfileLikelihood/data/Profiles_SMFromABM_clean.mat";
files.profile_from_data_file = "../ProfileLikelihood/data/Profiles_SMFromData_clean.mat";

%% admit samples
out = admitSampledABMParameters(files);


% %% interpolate at finer grid
% pars = {lattice_parameters.values};
% npars_abm = numel(pars);
% nx = 3; % probably want to choose odd so the 3 calculated values along each dimension are used
% par_grid = cell(npars_abm,1);
% for i = 1:npars_abm
%     par_grid{i} = linspace(pars{i}(1),pars{i}(end),nx);
% end
% 
% PG = cell(npars_abm,1);
% [PG{:}] = ndgrid(par_grid{:});
% 
% Vq = cell(npars_sm,2);
% for pi = 1:npars_sm
%     Vq{pi,1} = interpn(pars{:},squeeze(BS(pi,:,:,:,:,:,:,:,1)),PG{:});
%     Vq{pi,2} = interpn(pars{:},squeeze(BS(pi,:,:,:,:,:,:,:,2)),PG{:});
% end
% 
% abm_region_1_log = false(size(Vq{1})); % bounded by lambda, alpha, and K
% [lambda_start,lambda_end] = getProfileBounds(PLFromData.out{1}([1,end],:),threshold);
% % lambda_start = PLFromData.out{1}(1,1);
% % lambda_end = PLFromData.out{1}(1,end);
% for i = 1:size(PLFromData.out{1},2)
%     if PLFromData.out{1}(1,i)<lambda_start || PLFromData.out{1}(1,i)>lambda_end
%         continue;
%     end
%     temp = true(size(Vq{1}));
%     for pi = 1:3
%         temp = temp & Vq{pi,1} <= PLFromData.out{1}(pi,i) & Vq{pi,2} >= PLFromData.out{1}(pi,i);
%     end
%     abm_region_1_log = abm_region_1_log | temp;
% end

% %% store the parameter vectors by both restrictions
% I = cell(7,1);
% [I{:}] = ind2sub(size(abm_region_1_log),find(abm_region_1_log));
% 
% LP1 =  lattice_parameters;
% for i = 1:7
%     LP1(i).values = reshape(par_grid{i}(I{i}),[],1);
% end
% 
% %%
% if nx~=3
%     warning("abm_region_1_log does not correspond to the sampled 3^7 grid.")
% end
% save("ABMParamEstimates_FromProfile_WithK","LP1","abm_region_1_log")

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../../../SampleABMFns/")


