%DoCoolPlot2Stem
%Just plot the stem cells
%1/31/13 KAN

%cd Data/AA_0__BB_5/HypLocations
cd Data/AA_300__BB_2/HypLocations
XYZHyp = dlmread('XYZloc_t_1.txt');
XYZkt = dlmread('XYZkt_t_1.txt');
seedrate2 = 0;
symchange = 0.02;
DoMigration = 3;
ii = 1080;

load('voxelgrid_t_1.mat');

plane = find(XYZHyp(:,3) == 25);

XYZ = XYZHyp(:,:);
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

p=patch(isosurface(flipped==1,0));set(p,'facecolor','red' ,'edgecolor', 'none');
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

for jjjj =1:length(XYZ)  %Just plot stems %firstState:firstState+numStates
    %% Get positions of all cells of given type
    %stateGroup = find( (State(:,1) == jjjj) ); 
    %if length(stateGroup) > 0
    coords =  10*XYZ(jjjj,:); % n by 3 array of (x,y,z) coordinates
    
    %set colorbi
    colorbi = [1 1 0];
%     if jjjj == 1
%         colorbi = [1 0 0];
%     elseif jjjj == 2
%         colorbi = [0 0 1];
%     elseif jjjj == 3
%         colorbi = [0 1 0];
%     end
    
    %% Plot spheres at each position, color coding states
    for kkkk = 1:size(coords,1)
        h1 = surf(coords(kkkk,1)+spherex,coords(kkkk,2)+spherey,coords(kkkk,3)+spherez);

        set(h1,'FaceLighting','phong','AmbientStrength',0.6,...
            'FaceColor', colorbi, 'EdgeAlpha', 0.1); % costly to set?
        %alpha(h1, 0.1);
    end
    %end
end  

disp(strcat("plotting ",num2str(length(XYZkt))," kt cells"));
XYZkt2 = XYZkt(:,:);
for iiii = 1:length(XYZkt2)
    coords =  10*XYZkt2(iiii,:); % n by 3 array of (x,y,z) coordinates
    colorbi = [0 1 0];
    %% Plot spheres at each position, color coding states
    for kkkk = 1:size(coords,1)
        h1 = surf(coords(kkkk,1)+spherex,coords(kkkk,2)+spherey,coords(kkkk,3)+spherez);

        set(h1,'FaceLighting','phong','AmbientStrength',0.6,...
            'FaceColor', colorbi, 'EdgeAlpha', 0.1); % costly to set?
        %alpha(h1, 0.1);
    end
end

%% Adjust axes settings
%brighten(0.2);
light
%lighting phong;
box on;		
grid off;

axis equal; 
%set(gca,'position',[0 0 1 1]);

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
