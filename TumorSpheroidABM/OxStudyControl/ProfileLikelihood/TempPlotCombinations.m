% a quick script for plotting parameter combinations after doing profile
% likelihood of the surrogate model on exp data (main_data.m)

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
load("ProfileLikelihoods_DataRestricted.mat")

%% profiles
for pi = 1:3
    figure;
    plot(out{pi}(pi,:),out{pi}(4,:))
    xlabel(sprintf("Profiled par: %s",para_names{pi}));
    ylabel("\chi^2")
    set(gca,'FontSize',20)
    savefig(sprintf("Profile_%s",regexprep(para_names{pi},"\","")))
    print(sprintf("Profile_%s",regexprep(para_names{pi},"\","")),'-dpng')
end


%% combinations
% this is saved as SMParameterCombinations.fig
figure;
ci=0;
for i = 1:2
    for j = i+1:3
        ci = ci+1;
        subplot(2,3,r2c(2,3,[1,ci]))
        plot(out{i}(i,:),out{i}(j,:))
        xlabel(sprintf("Profiled par: %s",para_names{i}));
        ylabel(sprintf("Best fit for: %s",para_names{j}));
        subplot(2,3,r2c(2,3,[2,ci]))
        plot(out{j}(j,:),out{j}(i,:))
        xlabel(sprintf("Profiled par: %s",para_names{j}));
        ylabel(sprintf("Best fit for: %s",para_names{i}));
    end
end
