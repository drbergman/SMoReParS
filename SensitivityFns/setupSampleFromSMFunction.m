function out = setupSampleFromSMFunction(files, sm_functional, options)

arguments
    files struct % structure containing the files with the necessary info for this
    sm_functional function_handle
    options.par_names {mustBeText} = strings(1,numel(x))
    options.T dictionary = configureDictionary("string","function_handle")
    options.D dictionary = configureDictionary("string","prob.ProbabilityDistribution")
    options.nsamps double {mustBeInteger} = 100
    options.sum_fn function_handle = @mean
    options.use_profiles logical = true
end

vals = loadABMValues(files);
if options.use_profiles
    BS = loadBoundingSurfaces(files.profiles);
    out = @(x) sampleFromSMProfiles(x, BS, vals, sm_functional, options);

else
    OP = loadOptimalParameters(files);
    out = @(x) sampleFromSMBest(x, OP, vals, sm_functional, options);
end
end

function vals = loadABMValues(files)

if ~isfield(files,"data")
    load(files.profiles,"files")
end
load(files.data,"vals")

end


function BS = loadBoundingSurfaces(profile_file)

load(profile_file,"profiles","files");
n_sm_pars = size(profiles,1);
profiles = reshape(profiles,n_sm_pars,[]);
n_abm_vecs = size(profiles,2);


BS = zeros(n_sm_pars,n_abm_vecs,2);
threshold = chi2inv(0.95,n_sm_pars);
for i = 1:n_abm_vecs
    for j = 1:n_sm_pars
        [BS(j,i,1),BS(j,i,2)] = getProfileBounds(profiles{j,i}([j,end],:),threshold);
    end
end

load(files.data,"cohort_size")
BS = reshape(BS,[n_sm_pars,cohort_size,2]);

end

function P = loadOptimalParameters(files)

if ~isfield(files,"optimal_parameters")
    load(files.profiles,"files")
end
load(files.optimal_parameters,"P")

end