% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

cohort_name = "cohort_2303301105";
nsamps = 5;
phase_dependent_death = true;
% 
% if ~phase_dependent_death
%     load("data/OptimalParameters.mat")
% else
%     load("data/OptimalParameters_phase_dependent_death.mat")
% end
load("data/OptimalParameters_UnLinkedHill.mat");

Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"count_*","state2_prop_*");
load(sprintf("../../data/%s/output.mat",cohort_name),"ids","lattice_parameters");
load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");
tt = tracked.t';
I = randperm(numel(P)/size(P,1),nsamps);

chemo_dim = 8; % dimension along which chemo concentration varies; make sure this is still dim 8!!
if ~any(cohort_name==["cohort_2303231625","cohort_2303271138"]) && ~isequal(lattice_parameters(chemo_dim).path,["chemo_pars","concentration"])
    error("make sure the chemo concentration dim is still the 8th!")
end

fn_opts.phase_dependent_death = true; % does chemo death occur over entirety of each phase (true)? Or is it a one-time event during a phase and so it happens at a higher rate during shorter phases (false)?
fn_opts.link_phase_death_rates = false; % whether to link the two phases death rates
fn_opts.hill_activation = true; % if unlinked, use hill activation?


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


