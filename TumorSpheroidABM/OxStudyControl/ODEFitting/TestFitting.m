% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

save_figs = true;
cohort_name = "cohort_230124175743017";

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ODEFittingFns/")

par_names = ["\lambda","\alpha","K"];

nsamps = 10;
par_file = "data/OptimalParameters_Using_DependentStatesAndTimeSeries.mat";
data_file = sprintf("../../data/%s/summary_short.mat",cohort_name);

fn = @computeTimeSeries;
fn_opts.condition_on_previous = false;

f = testSMFitToABM(par_file,data_file,nsamps,fn,fn_opts,par_names);

if save_figs
    fig_names = ["SampleFitsOfSMToABM_Using_DependentStatesAndTimeSeries","BestSMParameterDistributions_Using_DependentStatesAndTimeSeries"];
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

