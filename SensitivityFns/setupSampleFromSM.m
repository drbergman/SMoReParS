function out = setupSampleFromSM(files, sm_functional, options)

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

if options.use_profiles


out = @(x) sampleFromSMProfiles(x, BS, vals, sm_functional, options);