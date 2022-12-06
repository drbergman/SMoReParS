clearvars;
addpath("~/Documents/MATLAB/myfunctions/")

% this will create figures for a particular set of conditions that create
% the story (or at least the exposition thereof) for what happened

%% cohort info and loading
cohort_names = ["cohort_221013152528022","cohort_221012070831969","cohort_221012222824126","cohort_349087890609125"];
for dfi = 1:numel(cohort_names)
    cohort_name = cohort_names(dfi);
    load(sprintf("../data/%s/%s.mat",cohort_name,cohort_name))

    %% plotting options
    options = struct();

    %% check for subcohorts and set them up
    cohort_size = size(ids);
    sample_ind = find(cohort_size==nsamps_per_condition);
    if length(sample_ind)>1
        error("need to specify the index over which samples vary")
    elseif ~isvector(cohort_size) && sample_ind < length(cohort_size)
        error("samples are expected to vary over the final dimension")
    end

    if ~exist("subcohort_indices","var")
        subcohort_indices = 1:ndims(ids);
        subcohort_indices = setdiff(subcohort_indices,sample_ind);
    end

    n_subcohorts = prod(size(ids,subcohort_indices));

    if exist("subcohort_title_names","var")
        options.subcohort_title_names = subcohort_title_names;
    end

    if exist("subcohort_folder_names","var")
        options.subcohort_folder_names = subcohort_folder_names;
    end

    options.reprint = false;
    %% create a story book for each

    ids = reshape(ids,n_subcohorts,nsamps_per_condition);
    for ci = 1:n_subcohorts
        storybookForSubcohort(cohort_name,ids(ci,:)',ci,options);
        close all
    end
end