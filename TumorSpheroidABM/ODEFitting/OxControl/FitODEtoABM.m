% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;

% p.lambda = 10000;
% p.alpha = 10;
% p.K = 1e3;
% p.delta = .1;
% p.g1_prop0 = 0.1;
p(1) = 24/19; % lambda
p(2) = 24/5; % alpha
p(3) = 1e3; % K

lb = [0;0;0];
ub = [Inf;Inf;1e4];

cohort_name = "cohort_2303301105";
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

%% load ABM data
C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size");
Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"ode_state*");
load(sprintf("../../data/%s/output.mat",cohort_name),"ids");
load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");
t_abm = tracked.t;
compare_every = 0.25 / 24;

%% restrict to only these time points (probably should make this an interpolation rather than this
tt = 0:compare_every:round(t_abm(end));
[~,tind] = min(abs(t_abm - tt),[],1);

Sum.ode_state_average = reshape(Sum.ode_state_average,numel(t_abm),2,[]);
Sum.ode_state_std = reshape(Sum.ode_state_std,numel(t_abm),2,[]);

%%
P = zeros([npars,prod(C.cohort_size)]);

%% find the best fit pars; since the purpose of this is to give ranges for the ODE parameters for profiling later, the obj fun is not so critical
Avg = Sum.ode_state_average(tind,:,:);
Std = Sum.ode_state_std(tind,:,:);
FF(1:size(P,2)) = parallel.FevalFuture;
for i = 1:size(P,2)
    data = Avg(:,:,i);
    data_std = Std(:,:,i);
    F = @(p) sum(((computeTimeSeries(p,tt)./data - 1)).^2,'all'); % relative difference
%     F = @(p) sum(((computeTimeSeries(p,tt) - data)).^2,'all'); % difference
%     F = @(p) sum(((computeTimeSeries(p,tt) - data)./data_std).^2,'all'); % z-scores
    FF(i) = parfeval(@(x) fmincon(F,x0,[],[],[],[],lb,ub,[],opts),1,x0);
%     P(:,i) = fmincon(F,x0,[],[],[],[],lb,ub,[],opts);
end

for i = 1:size(P,2)

    [idx,temp] = fetchNext(FF);
    P(:,idx) = temp;

    if mod(i,100)==0
        fprintf("Finished %d.\n",i)
    end

end

