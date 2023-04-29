function out = rawError(p,tt,D,fn,C,fn_opts,input_opts)

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

sim_data = fn(p,tt,C,fn_opts,D.A);

opts = defaultRawErrorOptions;
if nargin>=7 && ~isempty(input_opts)
    opts = overrideDefaultOptions(opts,input_opts);
end

if opts.assume_independent_time_series % no opts passed in or chose to do independent time series
    out = -sum(((sim_data - D.A)./D.S).^2,"all","omitnan");
    if ~opts.only_use_z_scores
        out = 0.5*out - 0.5*size(sim_data,2)*log(2*pi()) - 0.5*sum(log(D.S),"all");
    end
else
    if opts.only_use_z_scores
        out = sum(my_mvnpdf(sim_data,D.A,D.C,struct("only_use_z_scores",true)));
    else
        out = sum(my_mvnpdf(sim_data,D.A,D.C,struct("only_use_z_scores",false)));
    end
end

if opts.report_as_error
    out = -out;
end

end

function default_options = defaultRawErrorOptions

default_options.assume_independent_time_series = true; % assume that the time series produced by the SM are independent (crazy, right? but it's what I had initially assumed, so this is the default value)
default_options.only_use_z_scores = true; % whether to use the constant and SD terms from the LL for normal distributions, or (if false) just use the sum of z-scores
default_options.report_as_error = true; % whether to report the value as an error for optimization purposes

end