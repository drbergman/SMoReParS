% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../../../SurrogateModelFns/")

save_figs = true;
reprint = true;
file_types = ["fig","png"];
fig_names = ["SampleFitsOfSMToABM_New","BestSMParameterDistributions_New","RSSOfSMFitsToABM_New"];
resolution = '-r300';

files.optimal_parameters = "data/SMFitToABM_New.mat";
files.sm_fit_file = "data/SMFitToData_New.mat";
load(files.optimal_parameters,"cohort_name")

column_names = ["G1/S Counts","G2/M Counts"];

par_names = ["\lambda","\alpha","K"];

rss_color = [102, 51, 153]/255;
nsamps = 4;
files.data = sprintf("../../data/%s/summary_short.mat",cohort_name);

sm.fn = @(p,t,c,opts,data) computeTimeSeries(p,t,data,opts.condition_on_previous,[]);
sm.opts.condition_on_previous = false;
sm.custom_raw_error_fn = @customRawError;

opts.column_names = column_names;
opts.par_names = par_names;
opts.place_par_names = "xlabel";

opts.abm_vec_inds = [2;12];

opts.rss_on_log_scale = true;
opts.rss_orientation = "Horizontal";
opts.rss_normalization = "count";
opts.show_rss_smoothing = true;

opts.data_color = "black";
fit_color = lines(2);
fit_color = sqrt(prod(fit_color,1));
opts.fit_color = fit_color;

[f,I] = testSMFitToABM(files,nsamps,sm,opts);
f(1).Children = flip(f(1).Children);
f(1).Children = reshape(reshape(f(1).Children,2,4)',[],1);
close(f(1:2))

%% margins for distributions figure
res_factor = 1;
f(3);
f(3).Name = "RSSForSMoReGloS";
f(3).Units = "inches";
f(3).Position(3:4) = res_factor*[1,1];
ax = gca;
ax.Box = "off";
ax.Children(1).LineWidth = res_factor*0.5;
ax.FontSize = res_factor*8;
margin = struct("left",.4,"right",.01,"top",.01,"bottom",.31);
spacing = struct("horizontal",0.2,"vertical",0.15);
uniformAxisSpacing(ax',margin,spacing);


%% save figures
saveFigures(f(3),file_types=file_types,save_figs=save_figs,reprint=reprint,resolution=resolution)

%% reset path
rmpath("../../../SurrogateModelFns/")
