% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100
% the SM is the Jain 2017 one with the return to proliferation rate forced
% nonnegative

clearvars;

make_save = false;
save_fig_opts.save_figs = false;
save_fig_opts.reprint = true;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = "SMFitToDataWithArrestedCompartments";

addpath("../../../ODEFittingFns/")

addpath("~/Documents/MATLAB/myfunctions/")
fixed_pars = [];

% p = [lambda,alphaRP,theta,VT,V0,alphaP,kalpha,a,rho0,delta0,kdelta,b,alphaR]
% fixed_pars = "theta";
[p,lb,ub] = basePars(fixed_pars);

fixed_pars = sort(fixed_pars); % just to keep these in the same order for saving
for i = 1:numel(fixed_pars)
    save_fig_opts.fig_names = save_fig_opts.fig_names + "_" + fixed_pars(i);
end

fn = @computeTimeSeriesWithArrestedCompartments;
fn_opts = [];
weight_choice = "uniform";

% opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
opts = optimset('Display','off');

%% finish loading/setting up
npars = length(p);

switch weight_choice
    case "uniform"
        weights = ones(3,1);
    case "only 1"
        weights = [1;0;0];
    case "only 2"
        weights = [0;1;0];
    case "only 3"
        weights = [0;0;1];
    case "not 1"
        weights = ~[1;0;0];
    case "not 2"
        weights = ~[0;1;0];
    case "not 3"
        weights = ~[0;0;1];
    case "random"
        weights = rand(3,1);
        weights = weights/sum(weights);
end

% Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
% Data point for 5 hours taken out, since it is incommensurate.

load("data/ExperimentalData.mat","t","D","C");

%%
sm = struct("fn",fn,"opts",fn_opts);
F = @(p) arrayfun(@(i) getRawError(sm,p,t,D(i),C{i}),1:3)*weights;
tic
[P,fstar] = fmincon(F,p,[],[],[],[],lb,ub,[],opts);
toc
AIC = fstar + 2*(npars-numel(fixed_pars));

%%
if make_save
    file_name = "data/ODEFitToDataWithArrestedCompartments";
    for i = 1:numel(fixed_pars)
        file_name = file_name + "_" + fixed_pars(i);
    end
    save(file_name + ".mat","P","fstar","weights","fn_opts","lb","ub")
end
%%
f=figure;
tfull = linspace(0,3,100);
ax = gobjects(2,3);
for i = 1:3
    for j = 1:2
        ax(j,i) = subplot(2,3,r2c(2,3,[j,i])); hold on;
    end
end
for i = 1:3
    sim_data = fn(P,tfull,C{i},fn_opts);
    plot(ax(1,i),tfull,sim_data(:,1),"--","LineWidth",2,"DisplayName","Fit");
    plot(ax(2,i),tfull,sim_data(:,2),"--","LineWidth",2,"DisplayName","Fit");
    patch(ax(1,i),[t;flip(t)],[D(i).A(:,1)-D(i).S(:,1);flip(D(i).A(:,1)+D(i).S(:,1))],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
    plot(ax(1,i),t,D(i).A(:,1),"black","Marker","o","MarkerFaceColor","black","DisplayName","Data");
    patch(ax(2,i),[t;flip(t)],[D(i).A(:,2)-D(i).S(:,2);flip(D(i).A(:,2)+D(i).S(:,2))],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
    plot(ax(2,i),t,D(i).A(:,2),"black","Marker","o","MarkerFaceColor","black","DisplayName","Data");
    title(ax(1,i),sprintf("C = %3.2fuM",C{i}),"Interpreter","none")
end
xlabel(ax,"Time (d)")
ylabel(ax(1,:),"Scaled Cell Count")
ylabel(ax(2,:),"Prop in G2/M")
sgtitle(sprintf("AIC = %3.4f",AIC))
normalizeYLims(ax(2,:))

% legend([fit_curve;data_curve;data_patch],"Location","northwest","FontSize",22)

set(ax,"FontSize",20)

saveFigures(f,save_fig_opts)

%% remove paths
rmpath("../../../ODEFittingFns/")
