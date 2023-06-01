function raw_error = rawError(p,tt,D,fn,C,fn_opts,input_opts)

% This function will compute the raw error (raw meaning it may be weighted
% later based on number of data points, etc) of the surrogate model. It 
% will take in the time, parameters, and data, as well as a function to
% compute the data and optional function options.
%
% p: parameters for the surrogate model to use
% tt: time points at which to evaluate the ODE
% D: structure with
%   D.A = average of data
%   D.S = standard deviation of data
% fn: function that turn the parameters, time points, and options into data
% to compare against D and S
% C: conditions for this particular sim; these are just passed into fn; use
% this for conditions that vary wihtin one call to the objective function,
% e.g. dose amount
% fn_opts: any options that are needed by the function; use this for
% options that do NOT vary within one call to the objective function, e.g.
% options that specify the version of the SM

opts = defaultRawErrorOptions;
if nargin>6 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

if opts.resample
    if opts.t(end) > tt(end)
        warning("Sampled time values appear to go beyond observed time values: %3.2f > %3.2f.\n",opts.t(end),tt(end))
    end
    D.A = interp1(tt,D.A,reshape(opts.t,[],1));
    if opts.assume_independent_time_series
        D.S = interp1(tt,D.S,reshape(opts.t,[],1));
    else
        D.C = interp1(tt,D.C,reshape(opts.t,[],1));
    end
    tt = opts.t;
end

sim_data = fn(p,tt,C,fn_opts,D.A);

if opts.assume_independent_time_series % no opts passed in or chose to do independent time series
    raw_error = -sum(((sim_data - D.A)./D.S).^2,"all","omitnan");
    if ~opts.only_use_z_scores
        raw_error = 0.5*raw_error - 0.5*size(sim_data,2)*log(2*pi()) - 0.5*sum(log(D.S),"all");
    end
else
    if opts.only_use_z_scores
        raw_error = sum(my_mvnpdf(sim_data,D.A,D.C,struct("only_use_z_scores",true)));
    else
        raw_error = sum(my_mvnpdf(sim_data,D.A,D.C,struct("only_use_z_scores",false)));
    end
end

if opts.report_as_error
    raw_error = -raw_error;
end

end

function default_options = defaultRawErrorOptions

default_options.assume_independent_time_series = true; % assume that the time series produced by the SM are independent (crazy, right? but it's what I had initially assumed, so this is the default value)
default_options.only_use_z_scores = true; % whether to use the constant and SD terms from the LL for normal distributions, or (if false) just use the sum of z-scores
default_options.report_as_error = true; % whether to report the value as an error for optimization purposes

default_options.resample = false; % whether or not to resample output at different time points
default_options.t = []; % time points to resample at


end