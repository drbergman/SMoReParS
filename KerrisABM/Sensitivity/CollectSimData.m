% This script will collect the HPC-run sims that will be used to compute
% global sensitivity.

clearvars;

%% identify the data
f = dir("data/sims/*/Number*.txt");
load("data/MOATLHSSample.mat","points"); % points at which the MOAT samples were taken
load("../PostAnalysis/summary.mat","display_par_names")

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
    if length(temp) < 300
        error("This sim did not make it to the endpoint.")
    end
    N(i) = temp(end);
end

N = reshape(N,npoints,nfacs+1,[]); % final dimension is for samples (I first loop over all parameter values, then over samples)
N = mean(N,3);

%% compute mu*

ee = zeros(nfacs,npoints);
for i = 1:npoints
    base_val = N(i,end);
    ee(:,i) = (N(i,1:nfacs)-base_val)/0.5;
end

mu_star = mean(abs(ee),2);
sigma = std(ee,[],2);
[mu_star,order] = sort(mu_star,"descend");
sigma = sigma(order);
display_par_names = display_par_names(order);

save("data/GlobalSensitivityDirect","mu_star","sigma","display_par_names","npoints");



