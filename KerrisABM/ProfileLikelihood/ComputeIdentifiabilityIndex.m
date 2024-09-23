clearvars;

addpath("../../ProfileLikelihoodFns/")

% model_type = ["exponential","logistic","von_bertalanffy"];
model_type = ["von_bertalanffy"];

for i = 1:length(model_type)
    files.profiles = sprintf("data/ProfileLikelihoods_%s_resampled_clean.mat",model_type(i));
    indices = identifiabilityIndex(files);

    file_name = sprintf("data/IdentifiabilityIndex_%s_resampled_clean",model_type(i));
    if exist(file_name,"file") && ~overwrite_save_files
        continue
    end
    save(file_name,"files","indices")

end



rmpath("../../ProfileLikelihoodFns/")
