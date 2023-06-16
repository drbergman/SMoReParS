clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ProfileLikelihoodFns/")
addpath("../ODEFitting/")

profile_file = "data/Profiles_SMFromABM_New_clean.mat";

save_opts.save_figs = true;
save_opts.reprint = true;
save_opts.file_types = ["fig","png"];
save_opts.resolution = '-r1200';

opts = struct();
opts.abm_vec_inds = [2,12];
opts.LineWidth = 0.5;
opts.place_par_names = "xlabel";

sm_par_display_names = ["\lambda","\alpha","K"];
nsamps = 2;
%% color
fit_color = lines(2);
fit_color = sqrt(prod(fit_color,1));
opts.LineColor = fit_color;

%% test profile
[f,ax,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names,opts);

fac = 8;
for abm_par_ind = 1:2
    for sm_par_ind = 1:3
        g(abm_par_ind,sm_par_ind) = figure("Name",sprintf("SpecificProfile_%d_%s",opts.abm_vec_inds(abm_par_ind),regexprep(sm_par_display_names(sm_par_ind),"\","")));
        axt(abm_par_ind,sm_par_ind) = gca;
        copyobj(ax(abm_par_ind,sm_par_ind).Children,axt(abm_par_ind,sm_par_ind));
    
        for i = 1:numel(axt(abm_par_ind,sm_par_ind).Children)
            axt(abm_par_ind,sm_par_ind).Children(i).LineWidth = .5*fac;
        end
        g(abm_par_ind,sm_par_ind).Units = "inches";
        g(abm_par_ind,sm_par_ind).Position(3) = .25*fac;
        g(abm_par_ind,sm_par_ind).Position(4) = .25*fac;
    end
end
set(axt,"FontSize",8*fac,"XTick",[],"YTick",[]);
% set([axt(:).Children],"LineWidth",.5*fac)
% 
% %% set margins
% margin = struct("left",.1,"right",.05,"top",.02,"bottom",.15);
% spacing = struct("horizontal",.04,"vertical",.09);
% uniformAxisSpacing(ax,margin,spacing);

%% save figures
saveFigures(g,save_opts)

%% reset path
rmpath("../../../ProfileLikelihoodFns/")
rmpath("../ODEFitting/")
