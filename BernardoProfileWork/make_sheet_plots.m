%% PLOT PARAMETER FIT VS ABM PARAMETERS
clear all;
close all;

confidence_interval_results;

% SHEET PLOTS

% file path variables
save_file_path = ".\sheet_plots\";
plot_loc = ".\experimental_data\mu_P_vs_theta.fig";

% load mu_P vs theta plot
fig = openfig(plot_loc);
plt = findobj(fig, 'Type', 'Scatter');
MU = plt.XData;
THETA = plt.YData;
hold off;

% find points in the parameter sheets
p_div_steps = 100;
p_div_step_size = (abm_p_div(end) - abm_p_div(1)) / p_div_steps;

k = 1; % count the number of admissible parameters found

% for each value of mu_P, theta, see if they lie between the confidence
% sheets for some value of p_div, and div_lim

% loop through each value of mu_P, theta in the mu_P vs theta graph
for i = 1:length(MU)
    % loop through p_div
    for j = 1:p_div_steps  
        p_div = abm_p_div(1)+j*p_div_step_size;
        % loop through div_lim
        for m = 1:abm_div_lim(end)-abm_div_lim(1) 
            div_lim = abm_div_lim(1)+m;
            % find the interpolated values of mu_P and theta from the
            % sheet plots corresponding to p_div and div_lim
            mu_P_low = bilinearly_interpolate(abm_p_div, abm_div_lim, ...
                                              muP_fits_low, p_div, div_lim);
            mu_P_up = bilinearly_interpolate(abm_p_div, abm_div_lim, ...
                                             muP_fits_up, p_div, div_lim);
            theta_low = bilinearly_interpolate(abm_p_div, abm_div_lim, ...
                                               theta_fits_low, p_div, div_lim);
            theta_up = bilinearly_interpolate(abm_p_div, abm_div_lim, ...
                                              theta_fits_up, p_div, div_lim);
            % if the current value of mu_P and theta lie between the
            % interpolated values of mu_P and theta for this specific p_div
            % and div_lim, then add the points to the sheet
            if (mu_P_low < MU(i) && MU(i) < mu_P_up) && MU(i) < 0.41 && ...
                             (theta_low < THETA(i) && THETA(i) < theta_up) 
                X(k) = MU(i);
                Y(k) = THETA(i);
                P_DIV(k) = abm_p_div(1)+j*p_div_step_size;
                DIV_LIM(k) = abm_div_lim(1)+m;
                k = k+1;

                if abm_div_lim(1)+m == 9
                    abm_p_div(1)+j*p_div_step_size

                end
            end
        end
    end
end


% display figures
figure(1) % figure for p_div
p_name = "p_{div}";
scatter3(X, Y, P_DIV)
hold on;
plot3(MU, THETA, zeros(size(MU)))
xlabel("\mu_P");
ylabel("\theta");
zlabel(p_name);
savefig(figure(1), save_file_path + p_name + "sheet_plot.fig");
hold off;

figure(2) % figure for div_lim
p_name = "div_{lim}";
scatter3(X, Y, DIV_LIM)
hold on;
plot3(MU, THETA, zeros(size(MU)))
xlabel("\mu_P");
ylabel("\theta");
zlabel(p_name);
savefig(figure(2), save_file_path + p_name + "sheet_plot.fig");
hold off;
