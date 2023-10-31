function profiles = profileLikelihood(pbest,t,D,C,objfn_constants,profile_params,input_opts)

% profiles each parameter in pbest. pbest is the best fit at the point. the
% time series to compare against is given in tt, data, and data_std.
% para_ranges sets bounds for profiling. lb and ub set bounds for
% optimization while another parameter is being profiled. opts carries
% fmincon options. threshold sets a cutoff for when profiling can be
% stopped. save_all_pars can be passed in as true to save all parameter
% values, not just the profiled parameter value

% In each direction of pbest(i) (left and right), the ith parameter changes
% at 1% of pbest(i) for 10 steps. If it is still within para_ranges(i,:)
% and the value of the output is within F(pbest)+threshold, then continue
% in that direction until one of those conditions fails.

% will not assume that the pbest came from the same objective function as
% we are using here

opts = defaultProfileLikelihoodOptions;
if nargin > 6 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

m = numel(D); % number of conditions used

F = @(p) arrayfun(@(j) rawError(p,t,D(j),objfn_constants.fn,C{j},...
    objfn_constants.fn_opts,opts.raw_error_opts),1:m)*objfn_constants.weights;

%% make sure pbest is best
[pbest,val_at_center] = fmincon(F,pbest,profile_params.A,profile_params.b,[],[],profile_params.lb,profile_params.ub,[],profile_params.opts);
pbest = reshape(pbest,[],1); % make sure it is a column vector
npars = numel(pbest);

profiles = cell(npars,1);
for i = 1:npars
    if opts.save_all_pars
        profiles{i} = [pbest;val_at_center];
        par_ind = i;
    else
        profiles{i} = [pbest(i);val_at_center];
        par_ind = 1;
    end
    for dir = [-1,1] % move left and right along this parameter dimension
        [profiles{i},pbest,val_at_center] = profileInOneDirection(profiles{i},F,pbest,val_at_center,i,dir,profile_params,opts.save_all_pars);
    end
    [~,order] = sort(profiles{i}(par_ind,:),"ascend");
    profiles{i} = profiles{i}(:,order);
    if any(profiles{i}<0,"all")
        error("Some negative")
    end
end

end

function default_options = defaultProfileLikelihoodOptions

default_options.save_all_pars = true;
default_options.raw_error_opts = [];

end