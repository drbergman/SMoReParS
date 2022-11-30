function [best_fits, best_fit_vals, fixed_param_values, confidence_endpoints] = prof_like(param_ranges, param_names, ...
                                                                                           index, subdiv, times, cell_counts_mean, cell_counts_std, ...
                                                                                           para_names_file_save, ...
                                                                                           save_file_path, threshold)
%% This function uses the ODE model vbmodel to find the
%% profile likelihood of a specific parameter, namely the parameter given by
%% the input index in the param_ranges array.
%
%% Inputs:
%%   param_ranges          - Nx2 array of parameter ranges (N is the number of parameters)
%%   param_names           - strings representing parameter names in the same order as
%%                           param_ranges.
%%   index                 - index in the param_ranges array of the parameter that is 
%%                           fixed for the profile likelihood.
%%   subdiv                - number of subdivision of fixed parameter range for the
%%                           profile likelihood plots.
%%   cell_counts_mean      - array of means of the input data at each time point.
%%   cell_counts_std       - array of standard deviations of the input data at each time
%%                           point.
%%   param_names_file_save - parameter names to be written in the file save (here because
%%                           the parameter names may have latex typesetting which may not 
%%                           allowed in file names).
%%  save_file_path         - path to save the plots and other produced files.
%%  threshold              - confidence interval threshold for the input data.


% set up parameter names
fixed_p_name = param_names(index);
varied_p_names = param_names;
varied_p_names(index) = [];

% set up names for the saved plots
fixed_p_name_file_save = para_names_file_save(index);
varied_p_names_file_save = para_names_file_save;
varied_p_names_file_save(index) = [];

% specify step sizes for the optimization
steps = (1/subdiv) .* ( param_ranges(:,2) - param_ranges(:,1));
% initial guess of parameters for minimization
p_ic = mean(param_ranges,2);
% p_ic = 0.5 * (param_ranges(:,2) + param_ranges(:,1));
p_ic(index) = [];
% varying parameters lower and upper bounds
p_ranges = param_ranges;
p_ranges(index,:) = [];
% steps for the varying parameters
varied_p_steps = steps;
varied_p_steps(index) = [];

opts = optimset('Display','off');

for i = 0:subdiv
    % display progress
    if (mod(i, 100) == 0)
        disp("Fixing parameter " + fixed_p_name + " and currently in iteration " + i);
    end
    % find best fits
    %% for original vb model
    f = @(p, t) vbmodel_fix(p, cell_counts_mean(1), index, param_ranges(index,1) + i*steps(index), t)./cell_counts_std;
    %% for general vb model
    %f = @(p, t) gen_vbmodel_fix(p, index, param_ranges(index,1) + i*steps(index), t)./cell_counts_std;
    %%
    [p_fit, resnorm] = lsqcurvefit(f, p_ic, times, cell_counts_mean./cell_counts_std, p_ranges(:,1), p_ranges(:,2),opts);
    % store best fits
    best_fits(i+1,:) = p_fit;
    best_fit_vals(i+1) = resnorm;
    % set best fit as initial condition for next fit
    p_ic = p_fit;
    % store the values of the fixed parameter
    fixed_param_values(i+1) = param_ranges(index,1) + i*steps(index);
end

% find the confidence interval endpoints of the parameter
%res_threshold = (min(best_fit_vals) + threshold);
%confidence_endpoints = find_intesection(fixed_param_values, best_fit_vals, res_threshold);
confidence_endpoints = [0 0];

% create res norm plot
figure(1)
plot(fixed_param_values, best_fit_vals)
hold on;

% plot thershold line
%res_threshold_array = res_threshold * ones(length(fixed_param_values));
%plot(fixed_param_values, res_threshold_array, "red")

% add labels and plot title
title(fixed_p_name)
xlabel(fixed_p_name)
ylabel('Residual Norm')

% add legend
%lgd = legend("Residual Norm", "Confidence Threshold");

% save plot
name = append(save_file_path, fixed_p_name_file_save(1) , "_resnorm_plot.fig");
savefig(figure(1), name(1));
hold off

% create parameter relationship plots
for i = 1:length(varied_p_names)

    figure(i+1)
    % create plot
    scatter(fixed_param_values, best_fits(:,i))
    % add labels and plot title
    title(fixed_p_name)
    xlabel(fixed_p_name)
    ylabel(varied_p_names(i))
    
    % save plot
    name = append(save_file_path, fixed_p_name_file_save , "_vs_" , varied_p_names_file_save(i) , ".fig");
    savefig(figure(i+1), name(1));
end
end

