function out = admitSampledABMParameters(files,input_opts)

opts = defaultAdmitSampledABMParametersOptions;
if nargin > 1 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

%% load from files
load(files.abm_data_file,"cohort_size","vals","condition_dim");
P_ABM = load(files.profile_from_abm_file,"profiles");
P_data = load(files.profile_from_data_file,"profiles");

%% process data and profiles
if ~isempty(condition_dim) % if there was a condition varied across ABM sims, then ignore this one
    vals(condition_dim) = [];
end
npars_sm = size(P_ABM.profiles,1);
P_ABM.profiles = reshape(P_ABM.profiles,npars_sm,[]);
npoints = size(P_ABM.profiles,2);

assert(npars_sm==size(P_data.profiles{1},1)-1) % need the data profiles to save all parameter values, not just the profiled parameter

%% create bounding hypersurfaces
BS = zeros(npars_sm,npoints,2);
threshold = chi2inv(0.95,npars_sm);
for i = 1:npoints
    for j = 1:npars_sm
        [BS(j,i,1),BS(j,i,2)] = getProfileBounds(P_ABM.profiles{j,i}([j,end],:),threshold);
    end
end

%% admit parameter vectors
out = admitByMethod(vals,BS,P_data.profiles,opts);

end

function out = admitByMethod(abm_vals,BS,data_profile,opts)

sz = cellfun(@(x) size(x,1),abm_vals);
sz = reshape(sz,1,[]); % make sure it's a row vector
if numel(sz)==1
    sz = [sz,1]; % make sure the size at least 2 elements
end
npars_sm = numel(data_profile);
out = false(sz);

switch opts.admission_method
    case "single_best"
        % Use only the single_best SM pars from the data to select ABM parameters
        for i = 1:npars_sm % loop over SM parameter profiles
            [~,I] = min(data_profile{i}(end,:)); % select single_best SM parameter
            temp = true(1,prod(sz));
            for pi = 1:npars_sm
                temp = temp & BS(pi,:,1)<=data_profile{i}(pi,I) & BS(pi,:,2)>=data_profile{i}(pi,I);
            end
            out = out | reshape(temp,sz);
        end

    case "all_profiles"
        % loop over all profiles, using each point (intersecting across the
        % parameter value projections), and unioning across points and
        % parameters
        for i = 1:npars_sm % loop over SM parameter profiles
            for j = 1:size(data_profile{i},2) % loop over all points in SM space found for this parameter
                temp = true(1,prod(sz));
                for pi = 1:npars_sm % loop over SM parameter values at this point
                    temp = temp & BS(pi,:,1)<=data_profile{i}(pi,j) & BS(pi,:,2)>=data_profile{i}(pi,j);
                end
                out = out | reshape(temp,sz);
            end
        end

    case "all_profiles_resampled"
        % same as "all_profiles" but each profile will be resampled to
        % evenly space out the points and (probably) increase their density
        for i = 1:npars_sm % loop over SM parameter profiles
            profile = interp1(data_profile{i}(i,:),data_profile{i}',linspace(data_profile{i}(i,1),data_profile{i}(i,end),opts.nsamples))';
            for j = 1:opts.nsamples % loop over all sampled points in SM space for this parameter
                temp = true(1,prod(sz));
                for pi = 1:npars_sm % loop over SM parameter values at this point
                    temp = temp & BS(pi,:,1)<=profile(pi,j) & BS(pi,:,2)>=profile(pi,j);
                end
                out = out | reshape(temp,sz);
            end
        end

    otherwise
        error("""%s"" is not a method for admitting ABM parameters.\n",opts.admission_method)

end




end

function default_options = defaultAdmitSampledABMParametersOptions

default_options.admission_method = "single_best"; % how to admit sampled ABM parameter vectors
default_options.nsamples = 100; % if admission_method=="all_profiles_resampled", how many samples to use in profile

end
