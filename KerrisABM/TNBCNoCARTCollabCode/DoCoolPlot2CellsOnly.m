%DoCoolPlot2Stem
%Just plot the stem cells
%1/31/13 KAN


%cd Data/AA_1__BB_2/HypLocations


XYZall = dlmread('XYZloc_t_10.txt');
Stateall = dlmread('Stateh_KT1_Trial2_t10.txt');
seedrate2 = 0;
symchange = 0.02;
DoMigration = 3;
ii = 1080;

%plane = find(XYZall(:,3) == 25);

XYZ = XYZall(:,:);
Stemh = find(Stateall(:,1) == 1);
Progenh = find(Stateall(:,1) == 2);
%CCR5h = find(Stateall(:,3) == 1);
%KTh = find(Stateall(:,4) == 1);
%both = find(Stateall(:,1) == 1 & Stateall(:,3) == 1);

%Stateall(CCR5h,1) = 3;
%Stateall(KTh,1) = 4;
Stateall(Progenh,1) = 2;
Stateall(Stemh,1) = 1;
State = Stateall(:,:);

xmin = 0; % need better estimates for these
xmax = 500;  
ymin = 0;  
ymax = 500; 
firstState = 0;
numStates = 3;

figure

CoolPlotRadius = 1; % allows smaller than 1 so can see inside...
Coolcellradius = 0.6;
[spherex,spherey,spherez] = sphere;% (8);
spherex = spherex.*CoolPlotRadius.*Coolcellradius;
spherey = spherey.*CoolPlotRadius.*Coolcellradius;
spherez = spherez.*CoolPlotRadius.*Coolcellradius;
hold all;

for jjjj = 1:4 %Just plot stems %firstState:firstState+numStates
    %% Get positions of all cells of given type
    stateGroup = find( (State(:,1) == jjjj) ); 
    if length(stateGroup) > 0
    coords =  XYZ(stateGroup,:); % n by 3 array of (x,y,z) coordinates
    
    %set colorbi
    colorbi = [1 1 0];
    if jjjj == 1
        colorbi = [1 0 0];
    elseif jjjj == 2
        colorbi = [0 0 1];
    elseif jjjj == 3
        colorbi = [0 1 0];
    elseif jjjj == 4
        colorbi = [1 1 0];
    end
    
    %% Plot spheres at each position, color coding states

    for kkkk = 1:size(coords,1)
        h1 = surf(coords(kkkk,1)+spherex,coords(kkkk,2)+spherey,coords(kkkk,3)+spherez);

        set(h1,'FaceLighting','phong','AmbientStrength',0.6,...
            'FaceColor', colorbi, 'EdgeAlpha', 0.1); % costly to set?
        
        %++alpha(h1, 0.1);
    end
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
axis([0 50, 0 50, 0 50]); %set 2nd 4th and 6th vector elements to manually fix plot size 


%Save figure
    filenamep = ['seed',num2str(seedrate2), 'sym', num2str(100*symchange),'Mig', num2str(DoMigration),'Stems', '_ii',num2str(ii), '.jpeg'];
		%eval(['print -dpict ',filenamep]); % following uses low bit graphic: imwrite(xx,mm,filename)
		print( '-dbmp', filenamep); %I think he's using a mac file
