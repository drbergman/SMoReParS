function pcfTimeSeriesPlot(r,t,pcf)

pcf_limits = [min(pcf,[],'all','omitnan'),max(pcf,[],'all','omitnan')];

if pcf_limits(1)==pcf_limits(2)
    return;
end

one_color = 1*ones(1,3);
ngrid = 1001;
max_colors = 101; % max colors for each of positive and negative directions
positive_corr_n = round(ngrid*(pcf_limits(2)-1)/diff(pcf_limits));
negative_corr_n = ngrid-positive_corr_n;

% unique_pos_colors = floor((max_colors+1)*linspace(0,1,positive_corr_n))/(max_colors+1);
% unique_neg_colors = floor((max_colors+1)*linspace(0,1,negative_corr_n))/(max_colors+1);
unique_pos_color_green = floor((max_colors+1)*linspace(0,one_color(2),positive_corr_n))/(max_colors+1);
unique_pos_color_blue = floor((max_colors+1)*linspace(0,one_color(3),positive_corr_n))/(max_colors+1);
unique_neg_colors_red = floor((max_colors+1)*linspace(0,one_color(1),negative_corr_n))/(max_colors+1);
unique_neg_colors_green = floor((max_colors+1)*linspace(0,one_color(2),negative_corr_n))/(max_colors+1);

% positive_colors = [ones(positive_corr_n,1),flip(unique_pos_colors')*[1,1]];
% negative_colors = [unique_neg_colors'*[1,1],ones(negative_corr_n,1)];
positive_colors = [linspace(1,one_color(1),positive_corr_n)',flip(unique_pos_color_green'),flip(unique_pos_color_blue')];
negative_colors = [unique_neg_colors_red',unique_neg_colors_green',linspace(1,one_color(3),negative_corr_n)'];

cmap = [negative_colors;positive_colors];


% contourf(gca,r,t,pcf,ngrid,'EdgeColor','none')
imagesc(gca,r,t,pcf)
xlabel('r (lattice sites)'); ylabel('time (days)')
colorbar;
colormap(gca,cmap);
clim(gca,pcf_limits)
