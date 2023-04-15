% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

save_figs = false;
cohort_name = "cohort_230124175743017";

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ODEFittingFns/")

par_names = ["\lambda","\alpha","K"];

nsamps = 10;
par_file = "data/OptimalParameters_Using_optimizeSMParsFromABM.mat";
data_file = sprintf("../../data/%s/summary.mat",cohort_name);

fn = @computeTimeSeries;

f = testSMFitToABM(par_file,data_file,nsamps,fn,[],par_names);

if save_figs
    fig_names = ["SampleFitsOfSMToABM","BestSMParameterDistributions"];
    for i = 1:numel(f)
        if isempty(f(i).Name)
            f(i).Name = fig_names(i);
        end
        savefig(f(i),sprintf("figures/fig/%s",f(i).Name))
        print(f(i),sprintf("figures/png/%s",f(i).Name),"-dpng")
    end
end

