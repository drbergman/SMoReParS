% This script will create the MOAT sample points
clearvars;

nfacs = 4; % number of ABM parameters
Ns = 65; % number of points on each curve
Nr = 2; % number of curves for each parameter
omega_max = 8;

par_names = ["p_{div}";"s_{div}";"r_{mig}";"p_{lim}"];

s = linspace(-pi,pi,Ns+1);
s(end) = [];
sz = [nfacs,Ns,Nr,nfacs];
points = zeros(sz);
for fi = 1:nfacs
    omega = ones(1,nfacs);
    omega(fi) = omega_max;
    omega(setdiff(1:nfacs,fi)) = mod(0:(nfacs-2),4) + 1;
    for ri = 1:Nr
        phi = 2*pi*rand(1,nfacs);
        points(:,:,ri,fi) = efastG(s,omega',phi');
    end
end
points = reshape(points,nfacs,[]);

%% icdfs
D = makeMOATDistributions(par_names);
npoints = size(points,2);
for i = 1:npoints
    for j = 1:nfacs
        points(j,i) = icdf(D(par_names(j)),points(j,i));
    end
end

points(4,:) = round(points(4,:));
points = reshape(points,sz);

save("data/eFASTSample.mat","points")
% 
% 
% delta=1/npoints;
% coord = (0.5*delta):delta:1;
% points = zeros(npoints,nfacs+1,nfacs);
% for j=1:nfacs
%     points(:,nfacs+1,j) = coord(randperm(length(coord)));
% end
% 
% for i = 1:npoints
%     for j = 1:nfacs
%         points(i,j,:) = points(i,nfacs+1,:);
%         if points(i,nfacs+1,j)<=0.5
%             points(i,j,j) = points(i,nfacs+1,j) + 0.5;
%         else
%             points(i,j,j) = points(i,nfacs+1,j) - 0.5;
%         end
%     end
% end
% 
% points = reshape(points,[],nfacs);
% 
% D = makeMOATDistributions(par_names);
% 
% for i = 1:size(points,1)
%     for j = 1:nfacs
%         points(i,j) = icdf(D(par_names(j)),points(i,j));
%     end
% end
% 
% points(:,4) = round(points(:,4));
% 
% save("data/MOATLHSSample.mat","points")