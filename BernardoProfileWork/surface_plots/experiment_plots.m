clear all;
close all;

%% figure path
% mu_P theta gamma
param_name = "gamma";
cancer_type = "breast";
theta_val = "7_29501";
%open_fig_path = ".\" + param_name + "_abm_plot.fig";
open_fig_path = ".\" + param_name + "_abm_plot_bilinear.fig";
save_fig_path = ".\experiment_plots\theta_" + theta_val + "\";

%% param value - converted from mm^3 to number of cells
% 0.1265 7.29501
% 0.167  9.08503
% 0.2165 11.2596
% 0.2615 13.2295
% 0.3065 15.1956
% 0.3515 17.1592
% 0.3965 19.1211

conversion_factor = 2.5e5;
p_val = 7.29501 .* [1, 1, 1; 1, 1, 1; 1, 1, 1];
% convert from theta to gamma
%
p_val = 1 - p_val.^(-1);

% ABM PARAMETERS
abm_p_div = [0.05, 0.125, 0.245];
abm_div_lim = [8, 12, 15];

% confidence region vertical padding so the plot looks better
padding = 5;

% open figure
fig = openfig(open_fig_path, 'new');
hold on;
surf(abm_div_lim, abm_p_div, p_val,'FaceColor', 'y', 'EdgeColor', 'none')
legend('Confidence Lower Bound', 'Confidence Upper Bound', 'Profile Likelihood Value')
%savefig(fig, save_fig_path + cancer_type + "_" + param_name + "_abm_plot.fig");
savefig(fig, save_fig_path + cancer_type + "_" + param_name + "_bilinear_abm_plot_.fig");
