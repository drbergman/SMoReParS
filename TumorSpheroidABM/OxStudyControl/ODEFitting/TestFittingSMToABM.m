% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ODEFittingFns/")

save_fig_opts.save_figs = true;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = ["SampleFitsOfSMToABM_New","BestSMParameterDistributions_New","RSSOfSMFitsToABM_New"];
save_fig_opts.resolution = '-r1200';

files.optimal_parameters = "data/SMFitToABM_New.mat";
files.sm_fit_file = "data/SMFitToData_New.mat";
load(files.optimal_parameters,"cohort_name")

column_names = ["G1/S Counts","G2/M Counts"];

par_names = ["\lambda","\alpha","K"];

rss_color = [102, 51, 153]/255;
nsamps = 4;
files.data = sprintf("../../data/%s/summary_short.mat",cohort_name);

fn = @computeTimeSeries;
fn_opts.condition_on_previous = false;

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

f = testSMFitToABM(files,nsamps,fn,fn_opts,opts);

%% finish sample fitting figures
f(1).Units = "inches";
f(1).Position(3:4) = 2;
set(f(1).Children,"FontSize",8);
ypos = zeros(numel(f(1).Children),1);
for i = 1:numel(f(1).Children)
    ypos(i) = f(1).Children(i).Position(2);
    set(f(1).Children(i).Children,"LineWidth",0.5)
    yticks(f(1).Children(i),f(1).Children(i).YTick([1,end]))
end
% set x label
min_y = min(ypos);
I = ypos==min_y;
xlabel(f(1).Children(I),"Time (d)")

% move title up
max_y = max(ypos);
I = find(ypos==max_y);

% f(1).Children(I(1)).Title.Units = "normalized";
% f(1).Children(I(2)).Title.Units = "normalized";
% f(1).Children(I(1)).Title.Position(2) = 1.09;
% f(1).Children(I(2)).Title.Position(2) = 1.09;

%% color par histograms
colors = {[237,35,39]/255,[108,190,70]/255,[58,84,164]/255};
par_color = containers.Map(par_names,colors);
figure(f(2));
for i = 1:numel(f(2).Children)
    I = strcmp(f(2).Children(i).XLabel.String,par_names);
    if any(I)
        color = par_color(par_names(I));
        ax = f(2).Children(i);
        for j = 1:numel(ax.Children)
            if isa(ax.Children(j),'matlab.graphics.chart.primitive.Histogram')
                ax.Children(j).FaceColor = color;
                break
            end
        end
    end
end
drawnow
axt = subplot(2,2,4);
hold on;

%% yline in RSS
load(files.sm_fit_file,"fstar")
ax = f(3).Children;
if opts.rss_on_log_scale
    fstar = log(fstar);
end
if ~isfield(opts,"rss_orientation") || strcmpi(opts.rss_orientation,"Vertical")
    fstar_line = xline(ax,fstar,"LineWidth",2,"DisplayName","RSS of SM to Data");
else
    fstar_line = yline(ax,fstar,"LineWidth",2,"DisplayName","RSS of SM to Data");
end
fstar_line.Color = [1 69/255 0];
legend(ax,ax.Children(1),"Location","best")
ax.FontSize = 20;

for i = 1:numel(ax.Children)
    copyobj(ax.Children(i),axt); hold on
end
xlabel(axt,ax.XLabel.String)
ylabel(axt,ax.YLabel.String)
yticks(axt,ax.YTick)
yticklabels(axt,ax.YTickLabel)
for j = 1:numel(axt.Children)
    if isa(axt.Children(j),'matlab.graphics.chart.primitive.Histogram')
        axt.Children(j).FaceColor = rss_color;
    end
    if isa(axt.Children(j),'matlab.graphics.chart.primitive.Line')
        axt.Children(j).Color = rss_color;
    end
end

set(f(2).Children,"FontSize",8)
f(2).Units = "inches";
f(2).Position(3:4) = 2;


%% margins for distributions figure
margin = struct("left",.12,"right",.1,"top",.025,"bottom",.15);
spacing = struct("horizontal",0.2,"vertical",0.15);
ax = reshape(flip(f(2).Children),2,2);
uniformAxisSpacing(ax',margin,spacing);


%% save figures
saveFigures(f,save_fig_opts)

%% reset path
rmpath("../../../ODEFittingFns/")
