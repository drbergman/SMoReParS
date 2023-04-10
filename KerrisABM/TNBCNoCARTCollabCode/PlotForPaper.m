function[fighandle]=PlotForPaper(iteration, aa,ab,bb,het, tumor, vasc, cart, macro, hyp, ant)

%This function is basically the same as DoCoolPlotCombo,  but it allows us
%to control what is being plotted based on conditions. Essentially, it's
%like a toggle based on conditions. We can plot any of the following:
%tumor, vasc, cart, macrophages, hypoxia.

%BY: Tina G 6/24/2021
%[tumor cells, vasc, CART, macrophages, hypoxic cells, antigens?]
conditions= [tumor,vasc,cart,macro, hyp, ant]; %this array will determine what gets plotted: tumor, vasc, CART,
%macrophages, hypoxia and antigens
% IF het=1, do gradated. if het=0, do binary
iter = iteration;
AA = aa;
AB = ab;
BB = bb;
dirc = strcat('Data');
cd(dirc);
if(het==1)
    cd Gradated
elseif (het==0)
    cd Binary
end

dirc1=strcat('AA_',num2str(AA),'__AB_',num2str(AB),'/BB_',num2str(BB),'/OverTime');
cd(dirc1)
States=dlmread(strcat('State_CART',num2str(AA),'_Trial',num2str(BB),'_t',int2str(iter),'.txt'));
HypCell=dlmread(strcat('HypCell_t_',int2str(iter),'.txt'));
%XYZm=dlmread(strcat('XYZm_CART',num2str(AA),'_Trial',num2str(BB),'_t',int2str(iter),'.txt'));

try
    XYZh=dlmread(strcat('XYZh_t_',int2str(iter),'.txt'));
catch
    XYZh = [-1 -1 -1];
    disp('no hypoxic')
end
%XYZh=dlmread(strcat('XYZcart_t',int2str(iter),'.txt'));

%XYZm=dlmread(strcat('XYZm_CART20_Trial3_t',int2str(iter),'.txt'));
%cd Data/AA_0__BB_5/HypLocations
%iter=4;
XYZ = dlmread(strcat('XYZloc_t_',int2str(iter),'.txt'));

XYZhreal=XYZh*10;
nocart=0;
try
    XYZkt = dlmread(strcat('XYZcart_t_',int2str(iter),'.txt'));
catch
    XYZkt = [-1 -1 -1];
    nocart=1;
    disp('no CART');
end
seedrate2 = 0;
symchange = 0.02;
DoMigration = 3;
ii = 1080;

load(strcat('voxelgrid_t_',int2str(iter),'.mat'));

plane = find(XYZ(:,3) == 25);
firstState = 0;
numStates = 3;
f = figure('Name',strcat('AA_',num2str(AA),'__AB_',num2str(AB),'__BB_',num2str(BB),'__t_',num2str(iter)),'NumberTitle','off');

hold on
xlim([0 500])%Set lims so matlab does not graph over Mv
ylim([0 500])
zlim([0 500])
axis([0 500, 0 500, 0 500]); %set 2nd 4th and 6th vector elements to manually fix plot size
view(150,36);
%need to flip bc isosurface uses row,col,z while voxelgrid.Agent is XYZ
flipped = permute(voxelgrid.Agent,[2 1 3]);

%PLOT VASCULATURE
if conditions(2)==1
    p=patch(isosurface(flipped==1,0));
    set(p,'facecolor','red' ,'edgecolor', 'none');
end

daspect([1 1 1])
%rotate3d on;
view([60 30]);
camlight
lighting gouraud

CoolPlotRadius = 2; % allows smaller than 1 so can see inside...
Coolcellradius = 5;
[spherex,spherey,spherez] = sphere;
spherex = spherex.*CoolPlotRadius.*Coolcellradius;
spherey = spherey.*CoolPlotRadius.*Coolcellradius;
spherez = spherez.*CoolPlotRadius.*Coolcellradius;

CartPlotRadius = 1.5; % allows smaller than 1 so can see inside...
Cartcellradius = 5;
[cspherex,cspherey,cspherez] = sphere;
cspherex =cspherex.*CartPlotRadius.*Cartcellradius;
cspherey = cspherey.*CartPlotRadius.*Cartcellradius;
cspherez = cspherez.*CartPlotRadius.*Cartcellradius;
hold all;

XYZreal=10*XYZ;
%PLOT CANCER CELLS
%use states to split cancer cells into stems and progenitors
stems=find(States==1);
progs=find(States==2);
%Get coordinates of stem and progs separately; they will be plotted in the
%if statement below.
stemXYZ=XYZreal(stems,:);
progXYZ=XYZreal(progs,:);


%PLOT TUMOR CELLS

if conditions(1)==1
    if conditions(5)==0 %if we have no hypoxia do red stem and blue prog
for jjjj =1:size(XYZreal,1)
    disp(XYZ);
    %% Get positions of all cells of given type
    coords =XYZreal(jjjj,:); % n by 3 array of (x,y,z) coordinates
    %set colorbi
    
     %differentiate b etween stem and prog based on toggle   
    colorbi = [1 1 0];
    if States(jjjj) == 1
        colorbi = [1 0 0];
    elseif States(jjjj) == 2
        colorbi = [0 1 1];
   
    end
    %% Plot spheres at each position, color coding states
    for kkkk = 1:size(coords,1)
        h1 = surf(coords(kkkk,1)+spherex,coords(kkkk,2)+spherey,coords(kkkk,3)+spherez);

        set(h1,'FaceLighting','phong','AmbientStrength',0.6,...
            'FaceColor', colorbi, 'EdgeAlpha', 0.1); % costly to set?
    end
end

    elseif(conditions(5)==1) %if we have hypoxia, do blue normoxic and purple hypoxic
for jjjj =1:size(XYZreal,1)  
    %% Get positions of all cells of given type
    coords =XYZreal(jjjj,:); % n by 3 array of (x,y,z) coordinates
       colorbi = [1 1 0];
    if HypCell(jjjj) == 1
        colorbi=[0.4940, 0.1840, 0.5560];
    elseif HypCell(jjjj) == 0
        colorbi = [0 1 1];
   
    end
    %% Plot spheres at each position, color coding states
    for kkkk = 1:size(coords,1)
        h1 = surf(coords(kkkk,1)+spherex,coords(kkkk,2)+spherey,coords(kkkk,3)+spherez);

        set(h1,'FaceLighting','phong','AmbientStrength',0.6,...
            'FaceColor', colorbi, 'EdgeAlpha', 0.1); % costly to set?
    end
end       
        
    end
end

%PLOT CART
XYZkt2 = XYZkt(:,:);
XYZkt2 = XYZkt2*10;
if conditions(3)==1 && nocart==0
for jjjj =1:size(XYZkt2,1)  %Just plot stems %firstState:firstState+numStates
    %% Get positions of all cells of given type
    coords =XYZkt2(jjjj,:); % n by 3 array of (x,y,z) coordinates
    colorbi = [1 1 0];
    %% Plot spheres at each position, color coding states
    for kkkk = 1:size(coords,1)
        h1 = surf(coords(kkkk,1)+cspherex,coords(kkkk,2)+cspherey,coords(kkkk,3)+spherez);
        set(h1,'FaceLighting','phong','AmbientStrength',0.6,...
            'FaceColor', colorbi, 'EdgeAlpha', 0.1); % costly to set?
    end
end      
end

%PLOT MACROPHAGES

if conditions(4)==1
for jjjj =1:size(XYZm,1)  %Just plot stems %firstState:firstState+numStates
    %% Get positions of all cells of given type
    coords =XYZm(jjjj,:); % n by 3 array of (x,y,z) coordinates
    colorbi = [0.75, 0.75, 0];
    %% Plot spheres at each position, color coding states
    for kkkk = 1:size(coords,1)
        h1 = surf(coords(kkkk,1)+spherex,coords(kkkk,2)+spherey,coords(kkkk,3)+spherez);
        set(h1,'FaceLighting','phong','AmbientStrength',0.6,...
            'FaceColor', colorbi, 'EdgeAlpha', 0.1); % costly to set?
    end
end 
end


%PLOT   *JUST* HYPOXIA
if conditions(5)==1
 if(conditions(1)==0)
for jjjj =1:size(XYZhreal,1)  %Just plot stems %firstState:firstState+numStates
    %% Get positions of all cells of given type
    coords =XYZhreal(jjjj,:); % n by 3 array of (x,y,z) coordinates
       colorbi = [0.4940, 0.1840, 0.5560];
    %% Plot spheres at each position, color coding states
    for kkkk = 1:size(coords,1)
        h1 = surf(coords(kkkk,1)+spherex,coords(kkkk,2)+spherey,coords(kkkk,3)+spherez);

        set(h1,'FaceLighting','phong','AmbientStrength',0.6,...
            'FaceColor', colorbi, 'EdgeAlpha', 0.1); % costly to set?
    end
end         

 end
end

if conditions(6)==1
    AntigenHetPlot = dlmread(strcat('AntigenExp_t_',int2str(iter),'.txt'));
    AnitgenHetPlot = round(AntigenHetPlot,3);
    %AntigentHetPlot = [zeros(1,length(AntigenHetPlot));AntigenHetPlot;zeros(1,length(AntigenHetPlot))]; %gradate along green
    colormap(autumn(length(AntigenHetPlot)))
    %scatter3(XYZreal(:,1),XYZreal(:,2),XYZreal(:,3),40,AnitgenHetPlot,'filled','MarkerEdgeColor','k');
for jjjj =1:size(XYZreal,1)  %Just plot stems %firstState:firstState+numStates
    %% Get positions of all cells of given type
    coords =XYZreal(jjjj,:); % n by 3 array of (x,y,z) coordinates
       %colorbi = [0.4940, 0.1840, 0.5560];
       %colorbi=AntigenHetPlot;
       colormap(autumn(length(AntigenHetPlot)))
    %% Plot spheres at each position, color coding states
    for kkkk = 1:size(coords,1)
        h1 = surf(coords(kkkk,1)+spherex,coords(kkkk,2)+spherey,coords(kkkk,3)+spherez);
        set(h1,'FaceLighting','phong','AmbientStrength',0.6,...
             'EdgeAlpha', 0.1); % costly to set?
        colormap(autumn(length(AntigenHetPlot)))
        colorbar
    end
end    
    %colorbar
end
cd ..
%CHANGE NAME ACCORDINGLY
filename = ['PFP tum', 'AA',num2str(aa), '_AB', num2str(ab),'_BB', num2str(bb),'_t',num2str(iter), '.tiff'];
print(f,filename,'-dpng','-r400'); 
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
cd ../../../..
fighandle=f;
end