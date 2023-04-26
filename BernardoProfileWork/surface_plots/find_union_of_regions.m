clear all;
close all;

fig_path = ".\bw_regions\";

%% find union of regions
for i = 1:7
    img = imread(fig_path + "theta_breast_admissible_region_" + i + ".tiff");
    img = img == 255; % convert image pixel values to 0's and 1's
    % separate RGB components of the image
    if i == 1
        % create the black and white image that will illustrate the region
        img_bw = img;
    else
        % keep only the overlapping region between the old and new image;
        img_bw = ((img_bw == 1) & (img == 1));
    end
    imshow(img_bw)
end

img_bw = 1 - img_bw;
img_bw = 0.9*img_bw;
img_bw = 1 - img_bw;

imwrite(img_bw, fig_path + "union_admissible_region.tiff");
