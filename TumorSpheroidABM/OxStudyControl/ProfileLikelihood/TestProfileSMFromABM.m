clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

profile_file = "data/Profiles_SMFromABM_New_clean.mat";

save_opts.save_figs = true;
save_opts.reprint = false;
save_opts.file_types = ["fig","png"];
save_opts.fig_names = "SampleProfilesOfSMFromABM_New";
save_opts.resolution = '-r1200';

opts = struct();
opts.abm_vec_inds = [2,12];
opts.LineWidth = 0.5;
opts.place_par_names = "xlabel";

sm_par_display_names = ["\lambda","\alpha","K"];
nsamps = 4;
%% color
fit_color = lines(2);
fit_color = sqrt(prod(fit_color,1));
opts.LineColor = fit_color;

%% test profile
[f,ax,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);

if nsamps==2 % if just doing this for the axes on the box demo of how SP accepts parameters
    f.Name = "ExampleProfilesFromABMForShow_New";
end

f.Units = "inches";
f.Position(3) = 4;
f.Position(4) = 2;
set(f.Children,"FontSize",8);
set(ax(:,2:end),"YTick",[])
for i = 1:numel(ax)
    ax(i).XAxis.Label.FontWeight = "normal";
end

%% set margins
margin = struct("left",.1,"right",.05,"top",.02,"bottom",.15);
spacing = struct("horizontal",.08,"vertical",.09);
uniformAxisSpacing(ax,margin,spacing);

%% save figures
saveFigures(f,save_opts)

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../ODEFitting/")
