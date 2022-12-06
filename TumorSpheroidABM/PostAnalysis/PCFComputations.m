clearvars;
addpath("~/Lattice-PCF-4-MATLAB/src")


%% prepare inputs
cohort_names = ["cohort_349087890609125","cohort_221013152528022","cohort_221012222824126","cohort_221012070831969"];
options = defaultOptions();
options.stop_at_edge = true;
options.is_cross_pcf = true;
options.multiple_points_per_site = false; % whether multiple points can occupy the same lattice site
options.proportion_to_use = 0.03;
options.min_to_use = 4e2;

options.rr = linspace(0,26,27);

summarize = true;
load_options = ["default","tumor_to_active","high_antigen_to_immune","high_antigen_mut_to_immune","high_antigen_nonmut_to_immune","low_antigen_to_immune","low_antigen_mut_to_immune","low_antigen_nonmut_to_immune"];

for ci = 1:numel(cohort_names)
    for i = 1:numel(load_options)

        switch load_options(i)
            case "default"
                filename = sprintf("../data/%s/pcf.mat",cohort_names(ci));
            case "tumor_to_active"
                filename = sprintf("../data/%s/pcf_tumor_to_active.mat",cohort_names(ci));
            case "high_antigen_to_immune"
                filename = sprintf("../data/%s/pcf_high_antigen_to_immune.mat",cohort_names(ci));
            case "high_antigen_mut_to_immune"
                filename = sprintf("../data/%s/pcf_high_antigen_mut_to_immune.mat",cohort_names(ci));
            case "high_antigen_nonmut_to_immune"
                filename = sprintf("../data/%s/pcf_high_antigen_nonmut_to_immune.mat",cohort_names(ci));
            case "low_antigen_to_immune"
                filename = sprintf("../data/%s/pcf_low_antigen_to_immune.mat",cohort_names(ci));
            case "low_antigen_mut_to_immune"
                filename = sprintf("../data/%s/pcf_low_antigen_mut_to_immune.mat",cohort_names(ci));
            case "low_antigen_nonmut_to_immune"
                filename = sprintf("../data/%s/pcf_low_antigen_nonmut_to_immune.mat",cohort_names(ci));

            otherwise
                warning("no name given")
                continue;
        end

        if exist(filename,"file")
            continue;
        end

        out = pcfCohort(sprintf("../data/%s/%s.mat",cohort_names(ci),cohort_names(ci)),options,summarize,load_options(i));

        load(sprintf("../data/%s/%s.mat",cohort_names(ci),cohort_names(ci)),"tracked")

        out = reshape(out,size(tracked));
        clear tracked

        save(filename,"-v7.3")

    end
end
beep