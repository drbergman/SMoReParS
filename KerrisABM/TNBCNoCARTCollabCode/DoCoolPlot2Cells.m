%DoCoolPlot2Stem
%Just plot the stem cells
%1/31/13 KAN


xmin = 0; % need better estimates for these
xmax = 250;  
ymin = 0;  
ymax = 250; 
firstState = 0;
numStates = 2;

CoolPlotRadius = 1; % allows smaller than 1 so can see inside...
Coolcellradius = 0.6;
[spherex,spherey,spherez] = sphere;% (8);
spherex = spherex.*CoolPlotRadius.*Coolcellradius;
spherey = spherey.*CoolPlotRadius.*Coolcellradius;
spherez = spherez.*CoolPlotRadius.*Coolcellradius;
hold all;

for jjjj = 1 %Just plot stems %firstState:firstState+numStates
    %% Get positions of all cells of given type
    stateGroup = find(Agentmat (XYZ(:,1)>xmin) & (XYZ(:,1)<xmax) & ...
        (XYZ(:,2)>ymin) & (XYZ(:,2)<ymax) & (x == jjjj) ); 
    coords =  XYZ(stateGroup,:); % n by 3 array of (x,y,z) coordinates
    
    %set colorbi
    colorbi = [0 0 0];
    if jjjj == 1
        colorbi = [1 0 0];
    elseif jjjj == 2
        colorbi = [0 1 0];
    end
    
    %% Plot spheres at each position, color coding states
    for kkkk = 1:size(coords,1)
        h1 = surf(coords(kkkk,1)+spherex,coords(kkkk,2)+spherey,coords(kkkk,3)+spherez);

        set(h1,'FaceLighting','phong','AmbientStrength',0.6,...
            'FaceColor', colorbi, 'EdgeAlpha', 0.1); % costly to set?
        alpha(h1, 0.1);
    end
end  

%% Adjust axes settings
%brighten(0.2);
light
%lighting phong;
box on;		
grid off;

%axis equal; 
%set(gca,'position',[0 0 1 1]);

rotate3d on;		
%hold off;		
%axis equal; 
view(3);

%Save figure
    %filenamep = ['seed',num2str(seedrate2), 'sym', num2str(100*symchange),'Mig', num2str(DoMigration),'Stems', '_ii',num2str(ii), '.jpeg'];
		%eval(['print -dpict ',filenamep]); % following uses low bit graphic: imwrite(xx,mm,filename)
		%print( '-dbmp', filenamep); %I think he's using a mac file
