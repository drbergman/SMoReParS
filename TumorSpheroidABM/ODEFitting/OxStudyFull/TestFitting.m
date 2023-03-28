% a quick script to test the fitting of the ODE parameters to the ABM data
% they fit

clearvars;

addpath("~/Documents/MATLAB/myfunctions/")

cohort_name = "cohort_2303231625";
nsamps = 5;
load("data/OptimalParameters.mat")
phase_dependent_death = true;

Sum = load(sprintf("../../data/%s/summary.mat",cohort_name),"ode_state*");
load(sprintf("../../data/%s/output.mat",cohort_name),"ids","lattice_parameters");
load(sprintf("../../data/sims/%s/output_final.mat",ids(1)),"tracked");
tt = tracked.t';
I = randperm(numel(P)/size(P,1),nsamps);

sz = size(ids);
sz(end) = [];

chemo_dim = 8; % dimension along which chemo concentration varies; make sure this is still dim 8!!
if cohort_name~="cohort_2303231625"
    error("make sure the chemo concentration dim is still the 8th!")
end

Avg = Sum.ode_state_average;
Avg = permute(Avg,[1,2,chemo_dim+2,setdiff(3:ndims(Avg),chemo_dim+2)]); % move chemo concentration dim to right after time and phase
Avg = reshape(Avg,numel(tt),2,size(Avg,3),[]); 

Std = Sum.ode_state_std;
Std = permute(Std,[1,2,chemo_dim+2,setdiff(3:ndims(Std),chemo_dim+2)]); % move chemo concentration dim to right after time and phase
Std = reshape(Std,numel(tt),2,size(Std,3),[]);

Sum.ode_state_average = reshape(Sum.ode_state_average,length(tt),2,[]);
Sum.ode_state_std = reshape(Sum.ode_state_std,length(tt),2,[]);
P = reshape(P,size(P,1),[]);

%%
figure;
ax = gobjects(nsamps,3);
for i = 1:nsamps
    for j = 1:3
        ax(i,j) = subplot(nsamps,3,r2c(nsamps,3,[i,j]),"NextPlot","add");           
    end
end
for i = 1:nsamps

    avg = Avg(:,:,:,I(i));
    std = Std(:,:,:,I(i));
    for j = 1:3
        for k = 1:2
            patch(ax(i,j),[tt;flip(tt)],[avg(:,k,j)-std(:,k,j);flipud(avg(:,k,j)+std(:,k,j))],"black","FaceAlpha",0.2,"EdgeColor","none")
            plot(ax(i,j),tt,avg(:,k,j),"black")
        end
        out = computeTimeSeries(P(:,I(i)),tt,lattice_parameters(chemo_dim).values(j),phase_dependent_death);
        plot(ax(i,j),tt,out,"--","LineWidth",2)
    end

end

%%
% par_names = ["lambda","alpha","K","delta","G1 prop_0"];
par_names = ["lambda","alpha","K","d"];
figure;
for i = 1:4
    subplot(4,1,i)
    histogram(P(i,:))
    title(par_names(i))
end

rmpath("~/Documents/MATLAB/myfunctions/")
