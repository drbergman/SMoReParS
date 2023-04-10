%DoCoolPlot2Stem
%Just plot the stem cells
%1/31/13 KAN

iter = 1;

for iter = 1:10

%cd Data/AA_0__BB_5/HypLocations
cd Data/AA_30__BB_1/HypLocations
XYZ = dlmread(strcat('XYZloc_t_',int2str(iter),'.txt'));
try
    XYZkt = dlmread(strcat('XYZcart_t_',int2str(iter),'.txt'));
catch
    XYZkt = [-1 -1 -1];
end
seedrate2 = 0;
symchange = 0.02;
DoMigration = 3;
ii = 1080;

load(strcat('voxelgrid_t_',int2str(iter),'.mat'));

plane = find(XYZ(:,3) == 25);

% XYZ = gpuArray(XYZ(:,:));
 %Stemh = find(Stateall(:,1) == 1);
 %CCR5h = find(Stateall(:,3) == 1);
 %both = find(Stateall(:,1) == 1 & Stateall(:,3) == 1);

 %Stateall(CCR5h,1) = 3;
%Stateall(both,1) = 4;


% xmin = 0; % need better estimates for these
% xmax = 500;  
% ymin = 0;  
% ymax = 500; 
firstState = 0;
numStates = 3;
         
voxelgrid.Agent(2,:,:) = voxelgrid.Agent(1,:,:)|voxelgrid.Agent(2,:,:);
voxelgrid.Agent(:,2,:) = voxelgrid.Agent(:,1,:)|voxelgrid.Agent(:,1,:);         
voxelgrid.Agent(:,:,2) = voxelgrid.Agent(:,:,1)|voxelgrid.Agent(:,:,2);          
voxelgrid.Agent(end-1,:,:) = voxelgrid.Agent(end,:,:)|voxelgrid.Agent(end-1,:,:);          
voxelgrid.Agent(:,end-1,:) = voxelgrid.Agent(:,end,:)|voxelgrid.Agent(:,end-1,:);          
voxelgrid.Agent(:,:,end-1) = voxelgrid.Agent(:,:,end)|voxelgrid.Agent(:,:,end-1);

figure
hold on
        xlim([0 500])%Set lims so matlab does not graph over Mv
        ylim([0 500])
        zlim([0 500])

%need to flip bc isosurface uses row,col,z while voxelgrid.Agent is XYZ
flipped = permute(voxelgrid.Agent,[2 1 3]);
%flipped(:,1) = voxelgrid.Agent(:,2);
%flipped(:,2) = voxelgrid.Agent(:,1);
%temp=flipped(:,1);
%flipped(:,1)=flipped(:,2);
%flipped(:,2)=temp;
p=patch(isosurface(flipped==1,0));
set(p,'facecolor','red' ,'edgecolor', 'none');
daspect([1 1 1])
%rotate3d on;
view([60 30]);
camlight
lighting gouraud



CoolPlotRadius = 1; % allows smaller than 1 so can see inside...
Coolcellradius = 5;
[spherex,spherey,spherez] = sphere;% (8);
spherex = spherex.*CoolPlotRadius.*Coolcellradius;
spherey = spherey.*CoolPlotRadius.*Coolcellradius;
spherez = spherez.*CoolPlotRadius.*Coolcellradius;
hold all;

disp(strcat("plotting ",num2str(length(XYZ))," cells"));

XYZ=10*XYZ;
scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'filled','MarkerEdgeColor','k')

disp(strcat("plotting ",num2str(length(XYZkt))," kt cells"));
XYZkt2 = XYZkt(:,:);
XYZkt2 = XYZkt2*10;
scatter3(XYZkt2(:,1),XYZkt2(:,2),XYZkt2(:,3),'filled','MarkerEdgeColor','k')

%% Adjust axes settings
%brighten(0.2);
light
%lighting phong;
box on;		
grid off;

axis equal; 
%set(gca,'position',[0 0 1 1]);

xlabel('x');
ylabel('y');
rotate3d on;		
hold off;		
axis equal; 
axis([0 500, 0 500, 0 500]); %set 2nd 4th and 6th vector elements to manually fix plot size
view(150,36);

%Save figure
    filenamep = ['seed',num2str(seedrate2), 'sym', num2str(100*symchange),'Mig', num2str(DoMigration),'Hyp', '_ii',num2str(ii), '.jpeg'];
		%eval(['print -dpict ',filenamep]); % following uses low bit graphic: imwrite(xx,mm,filename)
		print( '-dbmp', filenamep); %I think he's using a mac file
cd ../../..

end
