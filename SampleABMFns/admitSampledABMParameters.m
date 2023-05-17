function out = admitSampledABMParameters(files,input_opts)

opts = defaultAdmitSampledABMParametersOptions;
if nargin > 1 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

%% load from files
load(files.abm_data_file,"cohort_size","vals","condition_dim");
P_ABM = load(files.profile_from_abm_file,"out");
P_data = load(files.profile_from_data_file,"out");

%% process data and profiles
if ~isempty(condition_dim) % if there was a condition varied across ABM sims, then ignore this one
    vals(condition_dim) = [];
end
npars_sm = size(P_ABM.out,1);
P_ABM.out = reshape(P_ABM.out,npars_sm,[]);
npoints = size(P_ABM.out,2);

%% create bounding hypersurfaces
BS = zeros(npars_sm,npoints,2);
threshold = chi2inv(0.95,npars_sm);
for i = 1:npoints
    for j = 1:npars_sm
        [BS(j,i,1),BS(j,i,2)] = getProfileBounds(P_ABM.out{j,i}([j,end],:),threshold);
    end
end
% BS = reshape(BS,[npars_sm,cohort_size,2]);

%% admit parameter vectors
out = admitByMethod(vals,BS,P_data.out,opts.admission_method);

end

function out = admitByMethod(abm_vals,BS,data_profile,admission_method)

switch admission_method
    case "best"
        % Use only the best SM pars from the data to select ABM parameters
        sz = cellfun(@numel,abm_vals);
        sz = reshape(sz,1,[]); % make sure it's a row vector
        if numel(sz)==1
            sz = [sz,1]; % make sure the size at least 2 elements
        end
        npars_sm = numel(data_profile);
        out = false(sz);
        for i = 1:npars_sm % loop over SM parameter profiles
            [~,I] = min(data_profile{i}(end,:)); % select best SM parameter
            temp = true(1,prod(sz));
            for pi = 1:npars_sm
                temp = temp & BS(pi,:,1)<=data_profile{i}(pi,I) & BS(pi,:,2)>=data_profile{i}(pi,I);
            end
            out = out | reshape(temp,sz);
        end
end




end

function default_options = defaultAdmitSampledABMParametersOptions

default_options.admission_method = "best"; % how to admit sampled ABM parameter vectors

end
