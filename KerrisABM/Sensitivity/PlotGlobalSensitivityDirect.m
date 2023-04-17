clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

filename = "GlobalSensitivityDirect";

line_width = 1;

load("data/GlobalSensitivityDirect.mat","mu_star","display_par_names","sigma","npoints")

c = categorical(display_par_names,display_par_names);

f=figureOnRight;
hold on;
bar(c,mu_star,"LineWidth",line_width)

std_error = sigma / sqrt(npoints);
er = errorbar(c,mu_star,zeros(size(std_error)),std_error,"CapSize",30);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
er.LineWidth = line_width;

ylabel("\mu^*")
title("Global Sensitivity: Direct Method")
set(gca,"FontSize",16)

if ~exist("figures/fig","dir")
    mkdir("figures/fig")
end
if ~exist("figures/png","dir")
    mkdir("figures/png")
end
savefig(sprintf("figures/fig/%s",filename))
print(sprintf("figures/png/%s",filename),"-dpng")

