function raw_error = rawError_FAV(p,tt,D,fn,C,fn_opts, ...
    assume_independent_time_series,only_use_z_scores,report_as_error, ...
    resample,resample_t)

% rawError using Function Argument Validation (FAV)

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

arguments
    % necessary parameters
    p (:,1) double % even if input as a row vector, it will be converted to a column
    tt (:,1) double % even if input as a row vector, it will be converted to a column
    D struct
    fn function_handle
    C % condition information can be any data type
    fn_opts struct

    % optional parameters
    assume_independent_time_series logical = true % assume that the time series produced by the SM are independent (crazy, right? but it's what I had initially assumed, so this is the default value)
    only_use_z_scores logical = true % whether to use the constant and SD terms from the LL for normal distributions, or (if false) just use the sum of z-scores
    report_as_error logical = true % whether to report the value as an error for optimization purposes
    resample logical = false % whether or not to resample output at different time points
    resample_t double = [] % time points to resample at
end

if resample
    if resample_t(end) > tt(end)
        warning("Sampled time values appear to go beyond observed time values: %3.2f > %3.2f.\n",resample_t(end),tt(end))
    end
    D.A = interp1(tt,D.A,reshape(resample_t,[],1));
    if assume_independent_time_series
        D.S = interp1(tt,D.S,reshape(resample_t,[],1));
    else
        D.C = interp1(tt,D.C,reshape(resample_t,[],1));
    end
    tt = resample_t;
end

sim_data = fn(p,tt,C,fn_opts,D.A);

if assume_independent_time_series % no opts passed in or chose to do independent time series
    raw_error = -sum(((sim_data - D.A)./D.S).^2,"all","omitnan");
    if ~only_use_z_scores
        raw_error = 0.5*raw_error - 0.5*size(sim_data,2)*log(2*pi()) - 0.5*sum(log(D.S),"all");
    end
else
    if only_use_z_scores
        raw_error = sum(my_mvnpdf(sim_data,D.A,D.C,struct("only_use_z_scores",true)));
    else
        raw_error = sum(my_mvnpdf(sim_data,D.A,D.C,struct("only_use_z_scores",false)));
    end
end

if report_as_error
    raw_error = -raw_error;
end

end
