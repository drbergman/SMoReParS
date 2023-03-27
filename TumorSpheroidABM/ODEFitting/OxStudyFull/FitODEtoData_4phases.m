% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

clearvars;

p(1) = 24/11; % lambda
p(2) = 24/8;
p(3) = 24/4; % alpha
p(4) = 24/1;
p(5) = 1e3; % K
p(6) = 0.09; % chemo-induced death probability per uM of drug

lb = [0;0;0;0;0;0];
ub = [Inf;Inf;Inf;Inf;1e4;1/7.55];

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

npars = length(p);


% Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
% Data point for 5 hours taken out, since it is incommensurate.
D = load("data/ExperimentalData.mat");

%%
P = zeros(npars,1);
[pstar,fstar] = fmincon(@(p) F(p,D.tt,D.doses,D.data,D.data_std),x0,[],[],[],[],lb,ub,[],opts);

%%
figure;
tfull = linspace(0,3,100);
ax = gobjects(3,1);
for i = 1:3
    ax(i) = subplot(1,3,i); hold on;
    yy = computeTimeSeries4Phases(pstar,tfull,D.doses(i));
    if i == 1
        fit_curve = plot(tfull,sum(yy,2),"--","LineWidth",2,"DisplayName","Fit");
    else
        plot(tfull,[sum(yy(:,1:2),2),sum(yy(:,3:4),2)],"--","LineWidth",2)
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

function out = F(p,tt,doses,data,data_std)

out = sum(((sum(computeTimeSeries4Phases(p,tt,doses(1)),2) - data{1})./data_std{1}).^2,'all');
for j = 2:3
    temp = computeTimeSeries4Phases(p,tt,doses(j));
    temp = [sum(temp(:,1:2),2),sum(temp(:,3:4),2)];
    out = out + sum(((temp-data{j})./data_std{j}).^2,'all');
end

end
