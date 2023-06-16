% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

clearvars;
file_name = "SMFitToData";
make_save = true;
save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = file_name;

opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

%%
[p,lb,ub] = basePars();
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
    [P,fstar] = fmincon(F,p,[],[],[],[],lb,ub,[],opts);

%%
figure;
tfull = linspace(0,3,100);
yy = computeTimeSeries(P,tfull);
fit_curve = plot(tfull,sum(yy,2),"--","LineWidth",2,"DisplayName","Fit");
hold on;
data_curve = plot(tt,data,"black","Marker","o","MarkerFaceColor","black","DisplayName","Data");
data_patch = patch([tt;flip(tt)],[data-data_std;flip(data+data_std)],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
xlabel("Time (d)")
ylabel("Scaled Cell Count")

legend([fit_curve;data_curve;data_patch],"Location","northwest","FontSize",22)

set(gca,"FontSize",20)
