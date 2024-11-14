clearvars;

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

profile_file = "data/Profiles_SMFromABM_LMS_2_clean.mat";

save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = "SampleProfilesOfSMFromABM_LMS_bounded";
save_fig_opts.resolution = "-r1200";

load("../ODEFitting/data/SMFitToData_LMS_bounded","fixed_pars","model_type")

opts = struct();
opts.place_par_names = "xlabel";
opts.LineWidth = 0.5;
% opts.abm_vec_inds = 1:3; % set this to fix a subset of ABM parameters to test profiles
sm_par_display_names = ["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"delta";"kdelta";"rho0"];

nsamps = 4;

%% set up parameter names
D = parameterDisplayNameDictionary(model_type);
[~,I] = setdiff(sm_par_display_names,fixed_pars);
sm_par_display_names = sm_par_display_names(sort(I));
for i = 1:numel(sm_par_display_names)
    sm_par_display_names(i) = D(sm_par_display_names(i));
end

%% make plots
[f,ax,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);

%% finish figure
f.Units = "inches";
f.Position(3) = 5;
f.Position(4) = 3;
set(ax,"FontSize",8)
yticks(ax(:,2:end),[])
for i = 1:numel(ax)
    ax(i).XAxis.Label.FontWeight = "normal";
end

%% margins for distributions figure
margin = struct("left",.06,"right",.02,"top",.0,"bottom",.12);
spacing = struct("horizontal",0.05,"vertical",0.08);
uniformAxisSpacing(ax,margin,spacing);

%% save figures
saveFigures(f,save_fig_opts)

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../ODEFitting/")

