clear all;
close all;

%% region RGB
rgb = [255, 127, 128];
gray_tol = 150;

%% save figure as jpeg figure
cancer = "breast";
parameters = ["mu_P", "theta"];
% 0.1265 7_29501 0.86292
% 0.167  9_08503 0.88993
% 0.2165 11_2596 0.91119
% 0.2615 13_2295 0.92441
% 0.3065 15_1956
% 0.3515 17_1592
% 0.3965 19_1211 0.94770

theta_val = "7_29501";
o_fig_path = ".\experiment_plots\theta_" + theta_val + "\";
png_s_fig_path = ".\experiment_plots\theta_" + theta_val + "\";

%%{
for i = 1:length(parameters)
    figure(i) = openfig(o_fig_path + cancer + "_" + parameters(i) + "_abm_plot.fig", "new");
    view(2)
    legend off
    grid off
    %title("Amissible Region for " + parameters(i));
    exportgraphics(figure(i), ...
                    png_s_fig_path + cancer + "_" + parameters(i) + "_abm_plot.tiff", ...
                    "Resolution", 600);
end
%}

%% find admissible region
for i = 1:length(parameters)
    img = imread(png_s_fig_path + cancer + "_" + parameters(i) + "_abm_plot.tiff");
    % separate RGB components of the image
    R = img(:,:,1);
    G = img(:,:,2);
    B = img(:,:,3);
    if i == 1
        % create the black and white image that will illustrate the region
        img_bw = (R == rgb(1) & G == rgb(2) & B == rgb(3)) | ...
                  (R < gray_tol & G < gray_tol & B < gray_tol);
    else
        % keep only the overlapping region between the old and new image
        temp_img = (R == rgb(1) & G == rgb(2) & B == rgb(3)) | ...
                    (R < gray_tol & G < gray_tol & B < gray_tol);
        img_bw = img_bw == 1 & temp_img == 1;
    end
end

% make the region of intersection and axis labels a dark gray
img_bw = 0.9*img_bw;
img_bw = 1-img_bw;

imwrite(img_bw, png_s_fig_path + "theta_" + theta_val + "_" + cancer + "_admissible_region.tiff");
%imwrite(img_bw, png_s_fig_path + "theta_" + theta_val + "_gamma_" + cancer + "_admissible_region.tiff");
