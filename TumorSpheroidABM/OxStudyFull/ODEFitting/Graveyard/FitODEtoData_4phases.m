% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

clearvars;
addpath("../../../ODEFittingFns/")
addpath("~/Documents/MATLAB/myfunctions/")

p = zeros(6,1);
p(1) = 24/11; % lambda
p(2) = 24/8;
p(3) = 24/4; % alpha
p(4) = 24/1;
p(5) = 1e3; % K
p(6) = 0.09; % chemo-induced death probability per uM of drug

p_unlinked = [0.1;0.1]; % chemo-induced death rates if they are unlinked
p_hill = [5;.1]; % [hill coefficient ; EC50] if unlinked

fn = @computeTimeSeries4Phases;
fn_opts.phase_dependent_death = false; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?
fn_opts.link_phase_death_rates = false; % whether to link the two phases death rates
fn_opts.hill_activation = true; % if unlinked, use hill activation?

ub = [Inf;Inf;Inf;Inf;1e4;1/7.55]; % don't let the probability ever exceed 1 (hence ub(5) = 1/7.55)

weight_choice = "random";

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%%
if ~fn_opts.link_phase_death_rates
    p = [p(1:5);p_unlinked];
    ub(7) = 1/7.55;
    if fn_opts.hill_activation
        p = [p;p_hill];
        ub(6:7) = 1; % in this case, p(6:7) is the limiting probability, not the per uM probability
        ub(8:9) = [Inf;Inf];
    end
end

npars = length(p);
lb = zeros(npars,1);

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
D = load("data/ExperimentalData.mat");

%%
P = zeros(npars,1);
sm = struct("fn",fn,"opts",fn_opts);
F = @(p) arrayfun(@(i) getRawError(sm,p,D.tt,[D.count(:,i),D.state2_prop(:,i)],[D.sigma_count(:,i),D.sigma_state2_prop(:,i)],D.doses(i)),1:3)*weights;
[pstar,fstar] = fmincon(F,p,[],[],[],[],lb,ub,[],opts);


% [pstar,fstar] = fmincon(@(p) F(p,D.tt,D.doses,D.data,D.data_std),p,[],[],[],[],lb,ub,[],opts);

%%
figure;
tfull = linspace(0,3,100);
ax = gobjects(2,3);
for i = 1:3
    for j = 1:2
        ax(j,i) = subplot(2,3,r2c(2,3,[j,i])); hold on;
    end
end
for i = 1:3
    sim_data = fn(pstar,tfull,D.doses(i),fn_opts);
    plot(ax(1,i),tfull,sim_data(:,1),"--","LineWidth",2,"DisplayName","Fit");
    plot(ax(2,i),tfull,sim_data(:,2),"--","LineWidth",2,"DisplayName","Fit");
    patch(ax(1,i),[D.tt;flip(D.tt)],[D.count(:,i)-D.sigma_count(:,i);flip(D.count(:,i)+D.sigma_count(:,i))],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
    plot(ax(1,i),D.tt,D.count(:,i),"black","Marker","o","MarkerFaceColor","black","DisplayName","Data");
    patch(ax(2,i),[D.tt;flip(D.tt)],[D.state2_prop(:,i)-D.sigma_state2_prop(:,i);flip(D.state2_prop(:,i)+D.sigma_state2_prop(:,i))],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
    plot(ax(2,i),D.tt,D.state2_prop(:,i),"black","Marker","o","MarkerFaceColor","black","DisplayName","Data");
    title(sprintf("C = %3.2fuM",D.doses(i)),"Interpreter","none")
end

xlabel(ax,"Time (d)")

% legend([fit_curve;data_curve;data_patch],"Location","northwest","FontSize",22)

set(ax,"FontSize",20)

% function out = F(p,tt,doses,data,data_std)
% 
% out = sum(((sum(computeTimeSeries4Phases(p,tt,doses(1)),2) - data{1})./data_std{1}).^2,'all');
% for j = 2:3
%     temp = computeTimeSeries4Phases(p,tt,doses(j));
%     temp = [sum(temp(:,1:2),2),sum(temp(:,3:4),2)];
%     out = out + sum(((temp-data{j})./data_std{j}).^2,'all');
% end
% 
% end

rmpath("../../../ODEFittingFns/")

