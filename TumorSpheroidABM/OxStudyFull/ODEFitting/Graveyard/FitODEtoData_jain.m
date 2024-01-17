% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

lambda = 5.96;
alphaRP = 0.99;
theta = 1;
VT = 548;
V0 = 55;
alphaP = 8.69;
kalpha = 21.67;
a = 0.5;
rho0 = 0.06;
delta0 = 0.96;
kdelta = 9.22;
b = 1;
alphaR = 3.0;
p = [lambda;alphaRP;theta;VT;V0;alphaP;kalpha;a;rho0;delta0;kdelta;b;alphaR];

fn = @computeTimeSeries_jain;
fn_opts = struct();
npars = length(p);

ub = [100;100;10;1e4;1e4;10;100;20;10;1e3;1e3;20;10];
lb = zeros(npars,1);
weight_choice = "uniform";

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%%

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

rmpath("../../../ODEFittingFns/")

