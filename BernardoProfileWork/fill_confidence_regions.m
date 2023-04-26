%% This script creates a gray background on the resnorm plots to illustrate
%% the confidence region.
clear all;
close all;

%% customizable variables
conf_interval = [9.98439	956.99];
param1 = "mu_p";
param2 = "theta";
open_fig_path = ".\plots\";
save_fig_path = ".\plots\05_8\with_filling\";

% creating the paths to open and save the figures with and without the
% confidence regions
fig_name = append(param1, "_vs_", param2, ".fig");
open_fig_path = open_fig_path + fig_name;
save_fig_path = save_fig_path + fig_name;

% confidence region color
gray_tone = 160 / 255;
color = [gray_tone gray_tone gray_tone];

% confidence region vertical padding so the plot looks better
padding = 5;

% open figure
fig = openfig(open_fig_path, 'new');
hold on;

% get y axis values to upper and lower bound to the shaded figure
fig = gcf;
ax_objs = fig.Children;
data_objs = ax_objs.Children;
y = data_objs(1).YData;
% find highest and lowest value of the y axis and add some room for the
% filled region
low_bd = min(y) - padding;
up_bd = max(y) + padding;
vert_bds = [low_bd up_bd up_bd low_bd];

% left and right bound for confidence region
hor_bds = [conf_interval(1) conf_interval(1) conf_interval(2) conf_interval(2)];

% create plot filling of the confidence intervals
p = fill(hor_bds, vert_bds, color,'FaceAlpha',0.3, ...
          'LineStyle','none');
% ensure the filling is bellow the plot
uistack(p,"bottom");

lgd = legend("Profile Likelihood Fits", "Confidence Region");
% save plot
savefig(fig, save_fig_path);
