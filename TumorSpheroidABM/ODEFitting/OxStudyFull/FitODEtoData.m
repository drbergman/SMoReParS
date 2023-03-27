% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

clearvars;

p(1) = 24/19; % lambda
p(2) = 24/5; % alpha
p(3) = 1e3; % K
p(4) = 0.1; % chemo-induced death rate per uM of drug
phase_dependent_death = true; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?

lb = [0;0;0;0];
ub = [Inf;Inf;1e4;2];

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

x0 = p;
npars = length(p);


% Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
% Data point for 5 hours taken out, since it is incommensurate.

D = load("data/ExperimentalData.mat");

%%
P = zeros(npars,1);

F = @(p) sum(((sum(computeTimeSeries(p,D.tt,D.doses(1),phase_dependent_death),2) - D.data{1})./D.data_std{1}).^2,'all') + ...
    sum(arrayfun(@(cci) sum(((computeTimeSeries(p,D.tt,D.doses(cci),phase_dependent_death) - D.data{cci})./D.data_std{cci}).^2,'all'),2:3));
[pstar,fstar] = fmincon(F,x0,[],[],[],[],lb,ub,[],opts);

%%
figure;
tfull = linspace(0,3,100);
ax = gobjects(3,1);
for i = 1:3
    ax(i) = subplot(1,3,i); hold on;
    yy = computeTimeSeries(pstar,tfull,D.doses(i),phase_dependent_death);
    if i == 1
        fit_curve = plot(tfull,sum(yy,2),"--","LineWidth",2,"DisplayName","Fit");
    else
        plot(tfull,yy,"--","LineWidth",2)
    end
    for j = 1:size(D.data{i},2)
        data_curve = plot(D.tt,D.data{i}(:,j),"black","Marker","o","MarkerFaceColor","black","DisplayName","Data");
        data_patch = patch([D.tt;flip(D.tt)],[D.data{i}(:,j)-D.data_std{i}(:,j);flip(D.data{i}(:,j)+D.data_std{i}(:,j))],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
    end
    title(sprintf("C = %3.2fuM",D.doses(i)),"Interpreter","none")
end
xlabel(ax,"Time (d)")
ylabel(ax,"Scaled Cell Count")

% legend([fit_curve;data_curve;data_patch],"Location","northwest","FontSize",22)

set(ax,"FontSize",20)
