
%% This script sets up and calls the profile likelihood method for
%% identifiability.

clearvars;
close all;

% folder to store plots and text files
save_file_path = "./s_shaped_data/";

% load data and compute the standard error
[time_points,cell_counts,cell_counts_mean,cell_counts_std,~] = data_mean_std();
time_points = time_points';

threshold = chi2inv(0.95,3); % compute threshold value for the parameter confidence intervals

% specify parameter ranges
para_ranges = [15e3	24e3;     % K
               1.56313	2.5;  % theta
               0.028	0.1];      % mu_p

% specify parameter names
para_names = {'K', '\theta', '\mu_p'};
para_names_file_save = {'K', 'theta', 'mu_p'};

N = size(para_ranges,1); % number of parameters
subdiv = 500; % parameter search step size

% file to store parameter confidence intervals
confidence_intervals_file_path = append(save_file_path, "confidence_intervals");
confidence_intervals_file = fopen(confidence_intervals_file_path,'w');

% store parameter ranges to a text file
para_ranges_file_path = append(save_file_path, "para_ranges");
para_ranges_file = fopen(para_ranges_file_path, 'w');

for i = 1:N
    %% do the profile likelihood search
    [best_fits_, fit_vals_, fixed_param_vals_, confidence_interval] = prof_like(para_ranges, ...
                                                            para_names, i, ...
                                                            subdiv, time_points, cell_counts_mean, cell_counts_std, ...
                                                            para_names_file_save, ...
                                                            save_file_path, ...
                                                            threshold);
  
    %% save values to files
    % get parameter name
    para_name = append(para_names_file_save(i), "");
    % store specified parameter range
    fprintf(para_ranges_file, '%s', para_name);
    fprintf(para_ranges_file, '\n');
    fprintf(para_ranges_file, '%g\t', para_ranges(i,:));
    fprintf(para_ranges_file, '\n');
    % store confidence interval
    fprintf(confidence_intervals_file, '%s', para_name);
    fprintf(confidence_intervals_file, '\n');
    fprintf(confidence_intervals_file, '%g\t', confidence_interval);
    fprintf(confidence_intervals_file, '\n');
end

% close files
fclose('all');

