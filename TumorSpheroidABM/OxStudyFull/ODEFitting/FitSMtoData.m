% fits the ODE parameters to the experimental data. data is scaled so the
% initial cell count is 100

clearvars;
file_name = "SMFitToData_LMS_bounded";
make_save = true;
save_fig_opts.save_figs = true;
save_fig_opts.reprint = false;
save_fig_opts.file_types = ["fig","png"];
save_fig_opts.fig_names = file_name;

addpath("../../../ODEFittingFns/")

addpath("~/Documents/MATLAB/myfunctions/")

model_type = "LogisticModelSimplified";
%         1     2   3   4      5      6    7   8       9  10  11
% p = [lambda,alpha,K,alphaR,alphaP,kalpha,a,delta0,kdelta,b,rho0]
fixed_pars = "rho0";
% fixed_pars = [];
fn = @computeTimeSeries;
[p,lb,ub,fn_opts.p_setup_fn,fixed_vals] = fixParameters(model_type,fixed_pars);

p = p.*exp(0.5*randn(size(p)));

weight_choice = "uniform";

optim_opts = optimset('Display','on','TolFun',1e-12,'TolX',1e-12,'MaxFunEval',1e4);

%%

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
F = @(p) arrayfun(@(i) rawError(p,t,D(i),fn,C{i},fn_opts),1:3)*weights;

[P,fstar] = fmincon(F,p,[],[],[],[],lb,ub,[],optim_opts);

%% save output
if make_save
    save("data/" + file_name,"P","fstar","weights","fn_opts","lb","ub","fixed_pars","fn","fn_opts","optim_opts","model_type","fixed_vals") %#ok<*UNRCH>
end
% %%
% f=figure;
% tfull = linspace(0,3,100);
% ax = gobjects(2,3);
% for i = 1:3
%     for j = 1:2
%         ax(j,i) = subplot(2,3,r2c(2,3,[j,i])); hold on;
%     end
% end
% for i = 1:3
%     sim_data = fn(P,tfull,C{i},fn_opts);
%     plot(ax(1,i),tfull,sim_data(:,1),"--","LineWidth",2,"DisplayName","Fit");
%     plot(ax(2,i),tfull,sim_data(:,2),"--","LineWidth",2,"DisplayName","Fit");
%     patch(ax(1,i),[t;flip(t)],[D(i).A(:,1)-D(i).S(:,1);flip(D(i).A(:,1)+D(i).S(:,1))],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
%     plot(ax(1,i),t,D(i).A(:,1),"black","Marker","o","MarkerFaceColor","black","DisplayName","Data");
%     patch(ax(2,i),[t;flip(t)],[D(i).A(:,2)-D(i).S(:,2);flip(D(i).A(:,2)+D(i).S(:,2))],"black","FaceAlpha",0.2,"EdgeColor","none","DisplayName","+/- SD");
%     plot(ax(2,i),t,D(i).A(:,2),"black","Marker","o","MarkerFaceColor","black","DisplayName","Data");
%     title(ax(1,i),sprintf("C = %3.2fuM",C{i}),"Interpreter","none")
% end
% xlabel(ax,"Time (d)")
% ylabel(ax(1,:),"Scaled Cell Count")
% ylabel(ax(2,:),"Prop in G2/M")
% normalizeYLims(ax(2,:))

% legend([fit_curve;data_curve;data_patch],"Location","northwest","FontSize",22)

% set(ax,"FontSize",20)
% 
% saveFigures(f,save_fig_opts)

%% remove paths
rmpath("../../../ODEFittingFns/")
