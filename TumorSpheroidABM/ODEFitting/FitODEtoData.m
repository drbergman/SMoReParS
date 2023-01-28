% clearvars;

% p.lambda = 10000;
% p.alpha = 10;
% p.K = 1e3;
% p.delta = .1;
% p.g1_prop0 = 0.1;
p(1) = 24/19;
p(2) = 24/5;
p(3) = 1e3;

lb = [0;0;0];
ub = [Inf;Inf;1e4];

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
%%
% fn = fieldnames(p);
% npars = numel(fn);
% x0 = zeros(npars,1);
% for i = 1:npars
%     x0(i) = p.(fn{i});
% end
x0 = p;
npars = length(p);


% Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
% Data point for 5 hours taken out, since it is incommensurate.

tt = [0      10     24     36     48     72    ]';        % hours
tt = tt/24;                                             % days

% Control data
data = [0.899  1.340  1.633  2.408  3.557  5.583]';   % millions of cells
data_std = [0.099  0.193  0.207  0.298  0.168  0.364]';   % millions of cells

%% scale data
factor = 100/data(1);
data_std = data_std * factor;
data = data * factor;
%%
P = zeros(npars,1);

%     F = @(p) sum(((computeTimeSeries(p,tt)./cell_count - 1)).^2,'all');
%     F = @(p) sum(((computeTimeSeries(p,tt) - data)).^2,'all');
    F = @(p) sum(((sum(computeTimeSeries(p,tt),2) - data)./data_std).^2,'all');
    [pstar,fstar] = fmincon(F,x0,[],[],[],[],lb,ub,[],opts);

%%
figure;
tfull = linspace(0,3,100);
yy = computeTimeSeries(pstar,tfull);
fit_curve = plot(tfull,sum(yy,2),"--","LineWidth",2,"DisplayName","Fit");
hold on;
data_curve = plot(tt,data,"black","Marker","o","MarkerFaceColor","black","DisplayName","Data");
data_patch = patch([tt;flip(tt)],[data-data_std;flip(data+data_std)],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
xlabel("Time (d)")
ylabel("Scaled Cell Count")

legend([fit_curve;data_curve;data_patch],"Location","northwest","FontSize",22)

set(gca,"FontSize",20)
