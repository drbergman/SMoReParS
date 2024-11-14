clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

profile_file = "data/Profiles_SMFromABM_New_clean.mat";

save_figs = false;
reprint = false;
file_types = ["fig","png"];
fig_names = "SampleProfilesOfSMFromABM_New_SMoReGloS_version";
resolution = '-r300';

opts = struct();
% opts.abm_vec_inds = [2,12];
opts.LineWidth = 0.5;
opts.place_par_names = "xlabel";

sm_par_display_names = ["\lambda","\alpha","K"];
nsamps = 4;
%% color
fit_color = lines(2);
fit_color = sqrt(prod(fit_color,1));
opts.LineColor = fit_color;
opts.abm_vec_inds = 1346;

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
    saveFigures(f,save_figs=save_figs,reprint=reprint,file_types=file_types,fig_names=fig_names,resolution=resolution)

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../ODEFitting/")
