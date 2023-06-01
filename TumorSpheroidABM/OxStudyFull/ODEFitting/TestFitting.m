% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

save_figs = true;
cohort_name = "cohort_2305311216";

addpath("~/Documents/MATLAB/myfunctions/")
addpath("../../../ODEFittingFns/")

par_names = ["\alpha_R";"\alpha_P";"k_\alpha";"\delta_0";"k_\delta";"\rho_0"];


nsamps = 10;
par_file = "data/OptimalParameters.mat";
data_file = sprintf("../../data/%s/summary.mat",cohort_name);

fn = @computeTimeSeries;

load("data/ODEFitToData.mat","fixed_pars")
D = parameterOrdering("LogisticModel");

[p,lb,ub] = basePars(fixed_pars);
fixed_inds = zeros(numel(fixed_pars),1);
for i = 1:numel(fixed_pars)
    fixed_inds(i) = D(fixed_pars(i));
end
lb(fixed_inds) = [];
ub(fixed_inds) = [];
fixed_vals = p(fixed_inds);
p(fixed_inds) = [];
fn_opts.p_setup_fn = @(p) this__p_setup_fn(p,fixed_inds,fixed_vals);

[f,I] = testSMFitToABM(par_file,data_file,nsamps,fn,fn_opts,par_names);

%% setup Avg and Std so that they go [time,death rate per uM,all other pars]
count = Sum.count_average;
count = permute(count,[1,chemo_dim+1,setdiff(2:ndims(count),chemo_dim+1)]); % move chemo concentration dim to right after time
count = reshape(count,numel(tt),size(count,2),[]); 

count_std = Sum.count_std;
count_std = permute(count_std,[1,chemo_dim+1,setdiff(2:ndims(count_std),chemo_dim+1)]); % move chemo concentration dim to right after time
count_std = reshape(count_std,numel(tt),size(count_std,2),[]);

state2_prop = Sum.state2_prop_mean;
state2_prop = permute(state2_prop,[1,chemo_dim+1,setdiff(2:ndims(state2_prop),chemo_dim+1)]); % move chemo concentration dim to right after time
state2_prop = reshape(state2_prop,numel(tt),size(state2_prop,2),[]); 

state2_prop_std = Sum.state2_prop_std;
state2_prop_std = permute(state2_prop_std,[1,chemo_dim+1,setdiff(2:ndims(state2_prop_std),chemo_dim+1)]); % move chemo concentration dim to right after time
state2_prop_std = reshape(state2_prop_std,numel(tt),size(state2_prop_std,2),[]);

P = reshape(P,size(P,1),[]);

%%
figure;
ax = gobjects(nsamps,3);
for i = 1:nsamps
    for j = 1:3
        ax(i,j) = subplot(nsamps,3,r2c(nsamps,3,[i,j]),"NextPlot","add");           
    end
end
figure;
ax_p = gobjects(nsamps,3);
for i = 1:nsamps
    for j = 1:3
        ax_p(i,j) = subplot(nsamps,3,r2c(nsamps,3,[i,j]),"NextPlot","add");           
    end
end
for i = 1:nsamps
    y = count(:,:,I(i));
    s = count_std(:,:,I(i));
    p = state2_prop(:,:,I(i));
    ps = state2_prop_std(:,:,I(i));
    for j = 1:3
        patch(ax(i,j),[tt;flip(tt)],[y(:,j)-s(:,j);flipud(y(:,j)+s(:,j))],"black","FaceAlpha",0.2,"EdgeColor","none")
        plot(ax(i,j),tt,y(:,j),"black")
        out = computeTimeSeries([P(1:5,I(i));3;P(6,I(i))],tt,lattice_parameters(chemo_dim).values(j),fn_opts);
        plot(ax(i,j),tt,out(:,1),"--","LineWidth",2)
        patch(ax_p(i,j),[tt;flip(tt)],[p(:,j)-ps(:,j);flipud(p(:,j)+ps(:,j))],"black","FaceAlpha",0.2,"EdgeColor","none")
        plot(ax_p(i,j),tt,p(:,j),"black")
        plot(ax_p(i,j),tt,out(:,2),"--","LineWidth",2)
    end

end

%%
% par_names = ["lambda","alpha","K","delta","G1 prop_0"];
par_names = ["lambda","alpha","K","d in G1/S","d in G2/M","Hill coefficient","EC50"];
figure;
for i = 1:size(P,1)
    subplot(size(P,1),1,i)
    histogram(P(i,:))
    title(par_names(i))
end

rmpath("../../../ODEFittingFns/")


function p = this__p_setup_fn(p_in,fixed_inds,fixed_vals)

p = zeros(11,1);
p(fixed_inds) = fixed_vals;
p(p==0) = p_in;

end