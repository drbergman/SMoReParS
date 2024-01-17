function raw_error = rawError(sm,p,tt,D,C,opts)

% This function will compute the raw error (raw meaning it may be weighted
% later based on number of data points, etc) of the surrogate model. It 
% will take in the time, parameters, and data, as well as a function to
% compute simulated data and optional function options.
%   * sm: structure of SM information that has the following fields:
%       * fn: function to simulate the SM
%       * opts: options to run with the SM. any options that are needed by the function; use this for
%           options that do NOT vary within one call to the objective function, e.g.
%           options that specify the version of the SM
%   * p: parameters for the surrogate model to use
%   * tt: time points at observations in D occur
%   * D: structure with data to compare against at times in tt. D.A is
%       required. One of D.S or D.C is required depending on
%       opts.assume_independent_time_series
%       * D.A = average of data
%       * D.S = standard deviation of data
%       * D.C = covariance matrix of data
%   * C: conditions for this particular sim; these are just passed into fn; use
%       this for conditions that vary wihtin one call to the objective function,
%       e.g. dose amount
%   * opts: options for this function. See defaultRawErrorOptions

%% process opts
opts = overrideDefaultOptions(defaultRawErrorOptions(),opts);

%% compute raw_error
if ~isempty(opts.resample_t)
    if opts.resample_t(end) > tt(end)
        warning("Sampled time values appear to go beyond observed time values: %3.2f > %3.2f.\n",opts.resample_t(end),tt(end))
    end
    D.A = interp1(tt,D.A,reshape(opts.resample_t,[],1));
    if opts.assume_independent_time_series
        D.S = interp1(tt,D.S,reshape(opts.resample_t,[],1));
    else
        D.C = interp1(tt,D.C,reshape(opts.resample_t,[],1));
    end
    tt = opts.resample_t;
end

sim_data = sm.fn(p,tt,C,sm.opts,D.A);

if opts.assume_independent_time_series
    raw_error = -sum(((sim_data - D.A)./D.S).^2,"all","omitnan");
    if ~opts.only_use_z_scores
        raw_error = 0.5*raw_error - 0.5*size(sim_data,2)*log(2*pi()) - 0.5*sum(log(D.S),"all");
    end
else
    raw_error = sum(my_mvnpdf(sim_data,D.A,D.C,struct("only_use_z_scores",opts.only_use_z_scores)));
end

if opts.only_use_z_scores
    raw_error = -raw_error;
end

end

function default_options = defaultRawErrorOptions()

default_options.assume_independent_time_series = true; % assume that the time series produced by the SM are independent (crazy, right? but it's what I had initially assumed, so this is the default value)
default_options.only_use_z_scores = true; % whether to use the constant and SD terms from the LL for normal distributions, or (if false) just use the sum of z-scores
default_options.report_as_error = true; % whether to report the value as an error for optimization purposes

default_options.resample_t = []; % time points to resample at (if not empty)

end
