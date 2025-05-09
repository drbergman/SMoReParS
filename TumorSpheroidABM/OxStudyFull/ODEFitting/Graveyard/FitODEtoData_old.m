% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

clearvars;

make_save = false;
save_fig_opts.save_figs = false;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = "SMFitToData";

addpath("../../../ODEFittingFns/")

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions

p = zeros(4,1);
p(1) = 24/19; % lambda
p(2) = 24/5; % alpha
p(3) = 1e3; % K
p(4) = 0.1; % chemo-induced death rate per uM of drug

p_unlinked = [1.9;1.9]; % chemo-induced death rates if they are unlinked
p_hill = [3;.1]; % [hill coefficient ; EC50] if unlinked

fn = @computeTimeSeries;
fn_opts.phase_dependent_death = true; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?
fn_opts.link_phase_death_rates = false; % whether to link the two phases death rates
fn_opts.hill_activation = true; % if unlinked, use hill activation?

ub = [Inf;Inf;1e4;2];

weight_choice = "uniform";

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%%
if ~fn_opts.link_phase_death_rates
    p = [p(1:3);p_unlinked];
    ub(5) = 2;
    if fn_opts.hill_activation
        p = [p;p_hill];
        ub(6:7) = [3;Inf];
    end
end

npars = length(p);
lb = zeros(npars,1);

if ~fn_opts.link_phase_death_rates
    if fn_opts.hill_activation
        lb(6) = 3;
    end
end

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

[P,fstar] = fmincon(F,p,[],[],[],[],lb,ub,[],opts);

%%
if make_save
    save("data/ODEFitToData.mat","P","fstar","weights","fn_opts","lb","ub")
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
normalizeYLims(ax(2,:))

% legend([fit_curve;data_curve;data_patch],"Location","northwest","FontSize",22)

set(ax,"FontSize",20)

saveFigures(f,save_fig_opts)

%% remove paths
rmpath("../../../ODEFittingFns/")
