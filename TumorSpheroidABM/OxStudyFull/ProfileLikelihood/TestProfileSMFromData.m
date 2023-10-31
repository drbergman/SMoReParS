clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = "SampleProfilesOfSMFromData_LMS_bounded_converted";
save_fig_opts.resolution = '-r1200';

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

opts.LineWidth = 0.5;
opts.place_par_names = "xlabel";
opts.parameter_layout = [3,3];
load("../ODEFitting/data/SMFitToData_LMS_bounded.mat","fixed_pars","model_type");
D = parameterDisplayNameDictionary(model_type);
sm_par_display_names = ["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"delta";"kdelta";"rho0"];

profile_file = "data/Profiles_SMFromData_LMS_bounded_clean_converted.mat";
nsamps = 1;

%% color
fit_color = lines(2);
fit_color = sqrt(prod(fit_color,1));
opts.LineColor = fit_color;

%% set up parameter names
[~,I] = setdiff(sm_par_display_names,fixed_pars);
sm_par_display_names = sm_par_display_names(sort(I));
for i = 1:numel(sm_par_display_names)
    sm_par_display_names(i) = D(sm_par_display_names(i));
end

%% test the profile
[f,ax,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);

%% prepare figure for paper
f.Units = "inches";
f.Position(3) = 2;
f.Position(4) = 2;
set(f.Children,"FontSize",8)

%% adjust margin
margin = struct("left",.1,"right",.08,"top",.02,"bottom",.18);
spacing = struct("horizontal",0.1,"vertical",0.19);
uniformAxisSpacing(ax,margin,spacing);

%% save the figures
saveFigures(f,save_fig_opts)

%% reset the path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../ODEFitting/")

