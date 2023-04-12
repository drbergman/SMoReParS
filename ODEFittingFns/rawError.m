function out = rawError(p,tt,D,fn,C,fn_opts)

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

sim_data = fn(p,tt,C,fn_opts);
out = sum(((sim_data - D.A)./D.S).^2,"all","omitnan");
