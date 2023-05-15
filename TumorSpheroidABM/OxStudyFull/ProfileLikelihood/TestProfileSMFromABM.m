clearvars;

save_figs = false;

addpath("../../../ProfileLikelihoodFns/")


sm_par_display_names = ["\lambda","\alpha","K","d_{G1/S}","d_{G2/M}","EC50"];
profile_file = "ProfileLikelihoods.mat";
nsamps = 5;
[f,I] = testProfileSMFromABM(profile_file,nsamps,sm_par_display_names);

if save_figs
    fig_names = "SampleProfilesOfSMFromABM";
    for i = 1:numel(f)
        if isempty(f(i).Name)
            f(i).Name = fig_names(i);
        end
        fig_folders = ["figures/fig","figures/png"];
        for j = 1:numel(fig_folders)
            if ~exist(fig_folders(j),"dir")
                mkdir(fig_folders(j))
            end
        end
        savefig(f(i),sprintf("figures/fig/%s",f(i).Name))
        print(f(i),sprintf("figures/png/%s",f(i).Name),"-dpng")
    end
end

rmpath("../../../ProfileLikelihoodFns/")

