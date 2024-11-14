clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions

model_type = "logistic";
% model_type = "von_bertalanffy";

suffix = "_very_large";

save_fig_opts.save_figs = false;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = sprintf("GlobalSensitivityMOATIndirect_%s%s",model_type,suffix);

line_width = 1;

load(sprintf("data/GlobalSensitivityMOATIndirect_%s%s.mat",model_type,suffix),"mu_star","display_par_names","sigma","npoints")

c = categorical(display_par_names,display_par_names);

f=figureOnRight;
hold on;
bar(c,mu_star,"LineWidth",line_width)

std_error = sigma / sqrt(npoints);
er = errorbar(c,mu_star,zeros(size(std_error)),std_error,"CapSize",30);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
er.LineWidth = line_width;

% set(gca,"YScale","log")
ylabel("\mu^*")
title("Global Sensitivity: Indirect Method" + this__title_fn(model_type))
set(gca,"FontSize",16)

saveFigures(f,save_fig_opts)

function out = this__title_fn(s)

s = regexprep(s,"_"," ");
s = strsplit(s);
for i = 1:length(s)
    temp = char(s(i));
    temp(1) = upper(temp(1));
    s(i) = string(temp);
end

out = " (" + strjoin(s) + ")";

end

