% This script will collect the HPC-run sims that will be used to compute
% global sensitivity.

clearvars;

%% identify the data
f = dir("data/sims/*/Number*.txt");
load("data/MOATLHSSample.mat","points"); % points at which the MOAT samples were taken
load("../PostAnalysis/data/summary.mat","display_par_names")

%% process points
nfacs = size(points,2);
points = reshape(points,[],nfacs+1,nfacs);
npoints = size(points,1);

%% process sim data
N = zeros(numel(f),1); % just get the final tumor size
for i = 1:numel(f)
    % Don't assume that dir returns the sims in the order I want
    I = strfind(f(i).folder,"_");
    num = str2double(f(i).folder(I+1:end));
    assert(num==i);
    % now I know this one is in the proper spot in the order

    temp = readmatrix(sprintf("%s/%s",f(i).folder,f(i).name));
    if length(temp) < 300 || temp(300)==0
        error("This sim did not make it to the endpoint.")
    end
    N(i) = temp(end);
end

N = reshape(N,npoints,nfacs+1,[]); % final dimension is for samples (I first loop over all parameter values, then over samples)
A = mean(N,3);
S = std(N,[],3);

%% compute mu*

ee = zeros(nfacs,npoints);
for i = 1:npoints
    base_val = A(i,end);
    ee(:,i) = (A(i,1:nfacs)-base_val)/0.5;
    decreased_ind = diag(squeeze(points(i,1:4,:))) - reshape(points(i,end,:),[],1) < 0;
    ee(decreased_ind,i) = -ee(decreased_ind,i);
end

mu_star = mean(abs(ee),2);
sigma = std(abs(ee),[],2);
[mu_star,order] = sort(mu_star,"descend");
sigma = sigma(order);
ordered_par_names = display_par_names(order);

% save("data/GlobalSensitivityMOATDirect","mu_star","sigma","ordered_par_names","npoints");



