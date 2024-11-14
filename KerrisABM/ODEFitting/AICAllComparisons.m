% A script to compute AIC scores for the SMs I'm using

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../SurrogateModelFns/")

save_figs = true;
reprint = true;
file_types = ["fig","png"];
resolution = '-r300';
% save_fig_opts.reprint_warning = false;
fig_names = ["AICComparisonBySM","ModelLikelihood"];


sm_model_palette = smModelPalette();

model_type = ["exponential","logistic","von_bertalanffy"];

data_file = "../PostAnalysis/data/summary.mat";

sm.fn = @computeTimeSeries;

resample = true;
resample_t = 15:15:75;


n_models = numel(model_type);
load(data_file,"t","D","C","cohort_size")
D = reshape(D,1,[]);
AIC = zeros([n_models,prod(cohort_size)]);
RSS = cell(n_models,1);
RSS_bar = cell(n_models,1);
for i = 1:n_models
    if model_type(i)=="von_bertalanffy"
        model_ode_fitting = load(sprintf("data/OptimalParameters_%s_resampled.mat",model_type(i)),"fstar","resample_t");
    else
        model_ode_fitting = load(sprintf("data/OptimalParameters_%s.mat",model_type(i)),"fstar","resample_t");
    end
    % model_ode_fitting = load(sprintf("data/OptimalParameters_%s.mat",model_type(i)),"fstar","resample_t");
    if isfield(model_ode_fitting,"resample_t")
        RSS{i} = model_ode_fitting.fstar;
        RSS_bar{i} = RSS{i} ./ length(model_ode_fitting.resample_t);
    else
        RSS{i} = model_ode_fitting.fstar;
        RSS_bar{i} = RSS{i} ./ 300;
    end
    AIC(i,:) = 2*i + RSS{i}(:);
end

%% make comparisons
comps = [2,1;3,1;3,2];

f = gobjects(size(comps,1),1);
for i = 1:3
    %% plot AIC by SM
    f(i)=figureOnRight("Name",sprintf("AIC_%s_%s",model_type(comps(i,1)),model_type(comps(i,2))));
    ax = gca;
    scatter(ax,AIC(comps(i,1),:),AIC(comps(i,2),:),'filled');
    xlabel(ax,this__title_fn(model_type(comps(i,1))));
    ylabel(ax,this__title_fn(model_type(comps(i,2))));
    xL = xlim;
    yL = ylim;
    axis(ax,"manual")
    M = max(xL(2),yL(2));
    line(ax,[0 M],[0 M])
    set(ax,"FontSize",16)
    title(ax,"AIC By SM")

    %% plot relative log-likelihood
    f(3+i) = figureOnRight("Name",sprintf("AIC_Waterfall_%s_%s",model_type(comps(i,1)),model_type(comps(i,2))));
    ax = gca;
    bar(sort(((AIC(comps(i,1),:)-AIC(comps(i,2),:))/2)))
    xlabel("ABM Parameter Vectors")
    ylabel("Relative Log-Likelihood")
    title("Relative Likelihood of SMs")
    % annotation("textarrow",.7*ones(1,2),[.85 .9],"String",sprintf("%s more likely",this__title_fn(model_types(comps(i,2)))))
    % annotation("textarrow",.7*ones(1,2),[.7 .65],"String",sprintf("%s more likely",this__title_fn(model_types(comps(i,1)))))
    set(ax,"FontSize",16)
end

%% all on single
f(end+1) = figureOnRight("Name","AllAICvVBScatter");
ax = gca;
hold on
vb_log = model_type=="von_bertalanffy";
rel_log_likelihood = (AIC(~vb_log,:) - AIC(vb_log,:));
rel_log_likelihood_pos_log = rel_log_likelihood > 0;
rel_log_likelihood_pos = rel_log_likelihood(rel_log_likelihood_pos_log);
log_rel_log_likelihood_pos = log10(rel_log_likelihood_pos);
rel_log_likelihood_neg = rel_log_likelihood(~rel_log_likelihood_pos_log);
log_rel_log_likelihood_neg = log10(-rel_log_likelihood_neg);
smallest_neg = min([log_rel_log_likelihood_pos;log_rel_log_likelihood_neg]);
smallest_neg_floor = floor(smallest_neg);
sn_diff = smallest_neg - smallest_neg_floor;
log_rel_log_likelihood_pos_relative_to_smallest = log_rel_log_likelihood_pos - smallest_neg_floor;
log_rel_log_likelihood_neg_relative_to_smallest = log_rel_log_likelihood_neg - smallest_neg_floor;
temp = zeros(size(rel_log_likelihood));
temp(rel_log_likelihood_pos_log) = log_rel_log_likelihood_pos_relative_to_smallest;
temp(~rel_log_likelihood_pos_log) = -log_rel_log_likelihood_neg_relative_to_smallest;
% temp2 = sign(rel_log_likelihood).*log10(abs(rel_log_likelihood));
% temp3 = rel_log_likelihood(:);
% temp3(temp3>0) = log10(temp3(temp3>0));
% temp3(temp3<0) = -log10(-temp3(temp3<0));
% temp3 = reshape(temp3,2,[]);
temp1 = temp(:,all(rel_log_likelihood>0,1)); % those where vB is best
temp2 = temp(:,rel_log_likelihood(1,:)>0 & rel_log_likelihood(2,:)<0); % those were log > vB > exp
temp3 = temp(:,all(rel_log_likelihood<0,1) & rel_log_likelihood(1,:)<rel_log_likelihood(2,:)); % those where exp > log > vB
temp4 = temp(:,all(rel_log_likelihood<0,1) & rel_log_likelihood(2,:)<rel_log_likelihood(1,:)); % those where log > exp > vB
% scatter(temp1(1,:),temp1(2,:),"filled","MarkerFaceColor",sm_model_palette("von_bertalanffy"))
% scatter(temp1(1,:),temp1(2,:),"filled","MarkerFaceColor","black")
% scatter(temp2(1,:),temp2(2,:),"filled","MarkerFaceColor",sm_model_palette("logistic"))

sz = 2;

% scatter(temp(1,[51,57,79]),temp(2,[51,57,79]),2*sz,"MarkerEdgeColor","black","MarkerFaceColor","none")
scatter(temp(1,57),temp(2,57),2*sz,"MarkerEdgeColor","black","MarkerFaceColor","none")

scatter(temp2(1,:),temp2(2,:),sz,"filled","MarkerEdgeColor",sm_model_palette("logistic"),"MarkerFaceColor",sm_model_palette("logistic"))
scatter(temp1(1,:),temp1(2,:),sz,"filled","MarkerEdgeColor",sm_model_palette("von_bertalanffy"),"MarkerFaceColor",sm_model_palette("von_bertalanffy"))
scatter(temp3(1,:),temp3(2,:),sz,"filled","MarkerEdgeColor",sm_model_palette("exponential"),"MarkerFaceColor",sm_model_palette("exponential"))
scatter(temp4(1,:),temp4(2,:),sz,"filled","MarkerEdgeColor",sm_model_palette("logistic"),"MarkerFaceColor",sm_model_palette("logistic"))

% scatter(temp2(1,:),temp2(2,:),"filled","MarkerFaceColor","black")
xL = xlim;
yL = ylim;
xlim(max(abs([xL,yL]))*[-1 1])
ylim(max(abs([xL,yL]))*[-1 1])
xL = xlim;
yL = ylim;
patch([xL(2), sn_diff, sn_diff, xL(2)],[yL(2) ,yL(2), sn_diff, sn_diff],sm_model_palette("von_bertalanffy"),"FaceAlpha",0.2,"EdgeColor","none")
patch([xL(1), -sn_diff, -sn_diff, xL(1)],[yL(2) ,yL(2), sn_diff, sn_diff],sm_model_palette("exponential"),"FaceAlpha",0.2,"EdgeColor","none")
patch([xL(2), sn_diff, sn_diff, xL(2)],[yL(1) ,yL(1), -sn_diff, -sn_diff],sm_model_palette("logistic"),"FaceAlpha",0.2,"EdgeColor","none")
patch([xL(1), -sn_diff, xL(1)],[yL(1), -sn_diff, -sn_diff],sm_model_palette("exponential"),"FaceAlpha",0.2,"EdgeColor","none")
patch([xL(1), -sn_diff, -sn_diff],[yL(1), -sn_diff, yL(1)],sm_model_palette("logistic"),"FaceAlpha",0.2,"EdgeColor","none")
line([xL(1), -sn_diff],[0 0],"LineStyle","--","Color","black")
line([sn_diff xL(2)],[0 0],"LineStyle","--","Color","black")
line([0 0],[yL(1), -sn_diff],"LineStyle","--","Color","black")
line([0 0],[sn_diff yL(2)],"LineStyle","--","Color","black")
ax.Children = flip(ax.Children);
xT = ceil(xL(1)):floor(xL(2));
xT(mod(xT,2)==1) = [];
xticks(xT);
xticklabels(10.^(abs(xT)+smallest_neg_floor).*sign(xT))
yT = ceil(yL(1)):floor(yL(2));
yT(mod(yT,2)==1) = [];
yticks(yT);
yticklabels(10.^(abs(yT)+smallest_neg_floor).*sign(yT))

%%
ax.FontSize = 8;
f(end).Units = "inches";
f(end).Position(3:4) = [2,1.5];

f = f(end);
%% save figures
saveFigures(f, save_figs=save_figs, reprint=reprint, file_types=file_types, resolution=resolution, fig_names=fig_names);

%% reset path
rmpath("../../SurrogateModelFns/")


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

