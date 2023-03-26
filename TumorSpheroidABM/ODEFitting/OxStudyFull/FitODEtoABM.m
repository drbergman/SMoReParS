% finds the best fit parameters for the ODE at each sampled ABM parameter
% vector. 

clearvars;

p(1) = 24/19; % lambda
p(2) = 24/5; % alpha
p(3) = 1e3; % K
p(4) = 0.1; % chemo-induced death rate per uM of drug

chemo_death_is_continuous = false; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?
lb = [0;0;0;0];
ub = [Inf;Inf;1e4;2]; % running this once showed that beyond 0.5, delta just decreases populations too fast

cohort_name = "cohort_2303231625";
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
C = load(sprintf("../../data/%s/output.mat",cohort_name),"cohort_size","lattice_parameters");
chemo_dim = 8; % dimension along which chemo concentration varies; make sure this is still dim 8!!
if cohort_name~="cohort_2303231625"
    error("make sure the chemo concentration dim is still the 8th!")
end

Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"ode_state*");
load(sprintf("../../data/%s/output.mat",cohort_name),"ids");
load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");
t_abm = tracked.t;

%% set chemo val
% chemo_val = zeros(size(ids));


%% setup Avg and Std so that they go [time,phase,death rate per uM,all other pars]
Avg = Sum.ode_state_average;
Avg = permute(Avg,[1,2,chemo_dim+2,setdiff(3:ndims(Avg),chemo_dim+2)]); % move chemo concentration dim to right after time and phase
Avg = reshape(Avg,numel(t_abm),2,size(Avg,3),[]); 

Std = Sum.ode_state_std;
Std = permute(Std,[1,2,chemo_dim+2,setdiff(3:ndims(Std),chemo_dim+2)]); % move chemo concentration dim to right after time and phase
Std = reshape(Std,numel(t_abm),2,size(Std,3),[]);

%%
P = zeros([npars,size(Avg,4)]); % one parameter vector that works for each concentration

%% find the best fit pars; since the purpose of this is to give ranges for the ODE parameters for profiling later, the obj fun is not so critical
FF(1:size(P,2)) = parallel.FevalFuture;
idx = cell(8,1);
sz = C.cohort_size(setdiff(1:length(C.cohort_size),chemo_dim));
nconc = C.cohort_size(chemo_dim);
for i = 1:size(P,2)
    data = Avg(:,:,:,i);
    data(data<=1) = 1; % when the data is too small (i.e. <=1), just compute the SM difference from 1
    % data_std = Std(:,:,:,i);
    F = @(p) sum(arrayfun(@(cci) sum((computeTimeSeriesChemo(p,t_abm,C.lattice_parameters(chemo_dim).values(cci),chemo_death_is_continuous)./data(:,:,cci) - 1).^2,'all'),1:nconc)); % relative difference
    FF(i) = parfeval(@(x) fmincon(F,x0,[],[],[],[],lb,ub,[],opts),1,x0);
end

for i = 1:size(P,2)

    [idx,temp] = fetchNext(FF);
    P(:,idx) = temp;

    if mod(i,100)==0
        fprintf("Finished %d.\n",i)
    end

end

