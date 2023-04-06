% in starting each profile, I used the best parameter value from a
% different objective function; this will go back to those best values and
% fix them

clearvars;

%% set up some parameters
cohort_name = "cohort_230124175743017";
compare_every = 6 / 24;
lb = [0;0;0];
ub = [Inf;Inf;1e4];
opts = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
npars = 3;

%% load best parameters
load("../../ODEFitting/OxControl/data/OptimalParameters_noapop.mat","P")

%% Load ABM output
C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size");
Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"ode_state*");
load(sprintf("../../data/%s/output.mat",cohort_name),"ids");
load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");

t_abm = tracked.t;
tt = 0:compare_every:round(t_abm(end));
[~,tind] = min(abs(t_abm - tt),[],1);

Sum.ode_state_average = reshape(Sum.ode_state_average,numel(t_abm),2,[]);
Sum.ode_state_std = reshape(Sum.ode_state_std,numel(t_abm),2,[]);

Avg = Sum.ode_state_average(tind,:,:);
Std = Sum.ode_state_std(tind,:,:);
%% load profiles
load("data/ProfileLikelihoods.mat","out")

%% set objective function
F = @(p,data,data_std) sum(((computeTimeSeries(p,tt) - data)./data_std).^2,'all'); % if the data has two phases

%% fix
sz = size(out);
figure;
hold on;
l = gobjects(2,1);
l(1) = plot(NaN,NaN,"Color","k","LineStyle","-","DisplayName","Original","LineWidth",2);
l(2) = plot(NaN,NaN,"Color","r","LineStyle",":","DisplayName","Updated","LineWidth",2);
show_plot = true;
for i = 1:numel(out)
    [ri,ci] = ind2sub(sz,i);
    data = Avg(:,:,ci);
    data_std = Std(:,:,ci);

    Aeq = zeros(1,npars);
    Aeq(ri) = 1;

    I = find(out{ri,ci}(1,:)==P(ri,ci));
    if numel(I)~=1
        error("should find just one")
    end
    if show_plot
        l(1).XData = out{ri,ci}(1,:);
        l(1).YData = out{ri,ci}(2,:);
    end
    [~,out{ri,ci}(2,I)] = fmincon(@(p) F(p,data,data_std),P(:,ci),[],[],Aeq,P(ri,ci),lb,ub,[],opts);
    if show_plot
        l(2).XData = out{ri,ci}(1,:);
        l(2).YData = out{ri,ci}(2,:);
        pause
    end
end