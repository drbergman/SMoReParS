% A script to compute AIC scores for the SMs I'm using

clearvars;

save_fig_opts.save_figs = true;
save_fig_opts.file_types = ["fig","png"];
% save_fig_opts.reprint_warning = false;
save_fig_opts.fig_names = ["AICComparisonBySM","ModelLikelihood"];

addpath("../../ODEFittingFns/")

model_types = ["von_bertalanffy","logistic"];

data_file = "../PostAnalysis/data/summary.mat";

sm.fn = @computeTimeSeries;

resample = true;
resample_t = 15:15:75;


n_models = numel(model_types);
load(data_file,"t","D","C","cohort_size")
D = reshape(D,1,[]);
AIC = zeros([n_models,prod(cohort_size)]);
for i = 1:n_models
    load(sprintf("data/OptimalParameters_%s.mat",model_types(i)),"P");
    P = reshape(P,size(P,1),[]);
    sm.opts.model_type = model_types(i);
    for j = 1:size(P,2)
        AIC(i,j) = getRawError(sm,P(:,j),t,D(j),C{1},resample = resample, resample_t = resample_t) + 2*size(P,1);
    end

end

%% plot AIC by SM
f=figureOnRight;
ax = gca;
scatter(ax,AIC(1,:),AIC(2,:),'filled');
xlabel(ax,this__title_fn(model_types(1)));
ylabel(ax,this__title_fn(model_types(2)));
xL = xlim;
yL = ylim;
axis(ax,"manual")
M = max(xL(2),yL(2));
line(ax,[0 M],[0 M])
set(ax,"FontSize",16)
title(ax,"AIC By SM")

%% plot relative log-likelihood
f(2) = figureOnRight;
ax = gca;
bar(sort(((AIC(1,:)-AIC(2,:))/2)))
xlabel("ABM Parameter Vectors")
ylabel("Relative Log-Likelihood")
title("Relative Likelihood of SMs")
annotation("textarrow",.7*ones(1,2),[.85 .9],"String",sprintf("%s more likely",this__title_fn(model_types(2))))
annotation("textarrow",.7*ones(1,2),[.7 .65],"String",sprintf("%s more likely",this__title_fn(model_types(1))))
set(ax,"FontSize",16)

%% save figures
saveFigures(f,save_fig_opts);

%% reset path
rmpath("../../ODEFittingFns/")


function out = this__title_fn(s)

s = regexprep(s,"_"," ");
s = strsplit(s);
for i = 1:length(s)
    temp = char(s(i));
    temp(1) = upper(temp(1));
    s(i) = string(temp);
end

out = strjoin(s);

end

