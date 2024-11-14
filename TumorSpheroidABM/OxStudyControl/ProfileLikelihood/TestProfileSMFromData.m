clearvars;

save_opts.save_figs = true;
save_opts.reprint = true;
save_opts.file_types = ["fig","png"];
save_opts.fig_names = "SampleProfilesOfSMFromData_New";
save_opts.resolution = '-r1200';

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

opts.LineWidth = 0.5;
opts.place_par_names = "xlabel";
% opts.abm_vec_inds = 1:3; % set this to fix a subset of ABM parameters to test profiles
sm_par_display_names = ["\lambda";"\alpha";"K"];

profile_file = "data/Profiles_SMFromData_New_clean.mat";
nsamps = 1;

%% color
fit_color = lines(2);
fit_color = sqrt(prod(fit_color,1));
opts.LineColor = fit_color;

%% test the profile
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);
ax = subplot(1,3,1);
ax(1).YLabel.String = "Best fit RSS";
ax(1).YLabel.FontWeight = "normal";

%% prepare figure for paper
f.Units = "inches";
f.Position(3) = 2;
f.Position(4) = 1;
ax = f.Children;
ax = flip(ax);
ax = ax';
set(ax,"FontSize",8)
for i = 1:numel(ax)
    ax(i).XAxis.Label.FontWeight = "normal";
end

%% adjust margin
margin = struct("left",.15,"right",.08,"top",.02,"bottom",.3);
spacing = struct("horizontal",0.09,"vertical",0.1);
uniformAxisSpacing(ax,margin,spacing);

%% save the figures
saveFigures(f,save_opts)

%% reset the path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../ODEFitting/")

