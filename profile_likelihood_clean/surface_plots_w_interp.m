close all;
clear all;

%% user defined variables
% import upper and lower sheet points
confidence_interval_results;
upper_bounds = muP_fits_up;
lower_bounds = muP_fits_low;
subdivisions = 1000; % how many interpolated points do we want between data points
para = "mu_P"; % parameter name for the surface plots

%% create interpolation points
n = length(abm_p_div);
grid_size = (n-1)*subdivisions+1;
p_div = zeros([1,grid_size]);
div_lim = zeros([1,grid_size]);

for i= 1:n-1
    delta_p_div = (abm_p_div(i+1) - abm_p_div(i))/subdivisions;
    delta_div_lim = (abm_div_lim(i+1) - abm_div_lim(i))/subdivisions;
    for j = 1:subdivisions
        index = (i-1)*subdivisions + j;
        p_div(index) = abm_p_div(i) + (j-1)*delta_p_div;
        div_lim(index) = abm_div_lim(i) + (j-1)*delta_div_lim;
    end
end
p_div(index+1) = abm_p_div(n);
div_lim(index+1) = abm_div_lim(n);


%% create interpolated surface points
F_low = zeros([(n-1)*subdivisions+1, (n-1)*subdivisions+1]);
F_up = zeros([(n-1)*subdivisions+1, (n-1)*subdivisions+1]);
% create upper and lower surface plots using bilinear interpolation
for i = 1:grid_size
    for j = 1:grid_size
        F_up(i,j) = bilinearly_interpolate(abm_p_div, abm_div_lim, upper_bounds, p_div(i), div_lim(j));
        F_low(i,j) = bilinearly_interpolate(abm_p_div, abm_div_lim, lower_bounds, p_div(i), div_lim(j));
    end
end

%% create surface plots
figure(1)
surf(div_lim, p_div, F_low, ...
     'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
%surf(div_lim, p_div, F_up, ...
%     'FaceColor', 'm', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%legend('Confidence Lower Bound', 'Confidence Upper Bound')
ylabel("ABM p_{div}");
xlabel("ABM div_{lim}");
zlabel(para)
hold off;
%savefig(figure(1), ".\surface_plots\"+para+"_abm_plot_with_edges.fig");
%savefig(figure(1), ".\surface_plots\"+para+"_abm_plot.fig");
savefig(figure(1), ".\surface_plots\lower_bound_"+para+"_abm_plot.fig");

