% This script will test the multi-dimensional profiles by sampling some of
% the surfaces found at a sampling of ABM parameter vectors, slicing at the
% mid value of each SM parameter

clearvars;

save_figure = false;

addpath("../../../ProfileLikelihoodFns/")
 
sm_par_display_names = ["\lambda","\alpha","K"];
profile_file = "data/MultiDimProfileLikelihoods.mat";
nsamps = 5;
opts.plot_type = "contourf";
opts.likelihood_scale = "linear";
opts.log_scale_pars = "K";
opts.transform_fn = @(x) exp((x-max(x,[],"all"))/abs(min(x-max(x,[],"all"),[],"all")));
[f,ax,I] = testMultiDimProfileSMFromABM(profile_file,nsamps,opts);


if save_figure
    file_name = "SampleMultiDimProfilesOfSMFromABM";
    fig_folders = ["fig","png"];
    for i = 1:numel(fig_folders)
        fig_folder_name = sprintf("figures/%s",fig_folders(i));
        if ~exist(fig_folder_name,"dir")
            mkdir(fig_folder_name)
        end
        file_path = sprintf("%s/%s",fig_folder_name,file_name);
        if fig_folders(i)=="fig"
            savefig(f,file_path)
        else
            print(f,file_path,sprintf("-d%s",fig_folders(i)))
        end
    end
end

rmpath("../../../ProfileLikelihoodFns/")

