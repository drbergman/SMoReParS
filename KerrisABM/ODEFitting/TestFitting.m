% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../ODEFittingFns/")

par_names = ["\alpha","\theta","\beta"];

nsamps = 36;
par_file = "data/OptimalParameters.mat";
data_file = "../PostAnalysis/summary.mat";

fn = @computeTimeSeries;

f = testSMFitToABM(par_file,data_file,nsamps,fn,[],par_names);

fig_names = ["SampleFitsOfSMToABM","BestSMParameterDistributions"];
for i = 1:numel(f)
    if isempty(f(i).Name)
        f(i).Name = fig_names(i);
    end
    savefig(f(i),sprintf("figures/fig/%s",f(i).Name))
    print(f(i),sprintf("figures/png/%s",f(i).Name),"-dpng")
end


