% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ODEFittingFns/")

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = ["SampleFitsOfSMToABM_LMS_bounded","BestSMParameterDistributions_LMS_bounded","RSSOfSMFitsToABM_LMS_bounded"];
save_fig_opts.resolution = '-r1200';

rss_color = [102, 51, 153]/255;
nsamps = 4;

files.optimal_parameters = "data/SMFitToABM_LMS_bounded.mat";

load(files.optimal_parameters,"cohort_name")

files.data = sprintf("../../data/%s/summary.mat",cohort_name);
files.sm_fit_file = "data/SMFitToData_LMS_bounded.mat";

load("data/SMFitToData_LMS_bounded.mat","fixed_pars","model_type");
D = parameterDisplayNameDictionary(model_type);
sm_par_display_names = ["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"delta";"kdelta";"rho0"];

% opts.column_names = {["Control","Count"],["Control","G2/M Fraction"];
%                 ["0.75\muM","Count"],["0.75\muM","G2/M Fraction"];
%                 ["7.55\muM","Count"],["7.55\muM","G2/M Fraction"]};
opts.column_names = {"Count","G2/M Fraction";
                "Count","G2/M Fraction";
                "Count","G2/M Fraction"};

opts.abm_vec_inds = [147,159,168,220];
opts.data_color = "black";
fit_color = lines(2);
fit_color = sqrt(prod(fit_color,1));
opts.fit_color = fit_color;
opts.place_par_names = "xlabel";
opts.rss_orientation = "Horizontal";

opts.rss_on_log_scale = true;
opts.rss_normalization = "count";
opts.show_rss_smoothing = true;

load("data/SMFitToData_LMS_bounded.mat","fn","fn_opts")

%% set up parameter names
[~,I] = setdiff(sm_par_display_names,fixed_pars);
sm_par_display_names = sm_par_display_names(sort(I));
for i = 1:numel(sm_par_display_names)
    sm_par_display_names(i) = D(sm_par_display_names(i));
end
opts.par_names = sm_par_display_names;

%% Test the fitting
[f,I] = testSMFitToABM(files,nsamps,fn,fn_opts,opts);

%% set up fit figure
f(1).Units = "inches";
f(1).Position(3) = 5;
f(1).Position(4) = 3;

ax = f(1).Children;
ax = flip(ax); % axes are naturally in reverse order of their creation
ax = reshape(ax,[6,nsamps])'; 
set(ax,"FontSize",8)
xlabel(ax(end,:),"Time (d)")

%% margins for fits figure
margin = struct("left",.075,"right",.04,"top",.05,"bottom",.105);
spacing = struct("horizontal",0.06,"vertical",0.06);
uniformAxisSpacing(ax,margin,spacing);


%% set up parameter histograms figure
f(2).Units = "inches";
f(2).Position(3) = 3;
f(2).Position(4) = 3;

ax = f(2).Children;
ax = flip(ax); % axes are naturally in reverse order of their creation
ax = reshape(ax,[3,3])'; 
set(ax,"FontSize",8)

set(ax,"YTick",[])

%% margins for parameter histograms
margin = struct("left",.015,"right",.02,"top",.01,"bottom",.12);
spacing = struct("horizontal",0.05,"vertical",0.14);
uniformAxisSpacing(ax',margin,spacing);

%% yline in RSS
load(files.sm_fit_file,"fstar")
ax = f(3).Children;
ax.Children(1).Color = rss_color;
ax.Children(2).FaceColor = rss_color;
if opts.rss_on_log_scale
    fstar = log(fstar);
end
if ~isfield(opts,"rss_orientation") || strcmpi(opts.rss_orientation,"Vertical")
    fstar_line = xline(ax,fstar,"LineWidth",2,"DisplayName","RSS of SM to Data");
else
    fstar_line = yline(ax,fstar,"LineWidth",2,"DisplayName","RSS of SM to Data");
end
fstar_line.Color = [1 69/255 0];
ax.FontSize = 8;
f(3).Units = "inches";
f(3).Position(3) = 2;
f(3).Position(4) = 3;

%% get xlim and labels right
if opts.rss_orientation=="Vertical"
    xL = xlim(ax);
    xL(1) = 0;
    xlim(ax,xL)
    xx=xticks(ax);
    xx = [0,xx];
    xticks(ax,xx)
    xticklabels(ax,['10^{0}',xticklabels(ax)])
else
    yL = ylim(ax);
    yL(1) = 0;
    ylim(ax,yL)
    yy=yticks(ax);
    yy = [0,yy];
    yticks(ax,yy)
    yticklabels(ax,['10^{0}',yticklabels(ax)])
end

%% margins for RSS figure
margin = struct("left",.18,"right",.05,"top",.01,"bottom",.11);
spacing = struct("horizontal",0.05,"vertical",0.14);
uniformAxisSpacing(ax',margin,spacing);

%% save the figures
saveFigures(f,save_fig_opts)

%% reset path
rmpath("../../../ODEFittingFns/")


