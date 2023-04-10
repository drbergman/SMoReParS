function[fighandle]=PlotEverything(iteration, aa,ab,bb,het, tumor, vasc, cart, macro, hyp, ant)
%INPUTS ARE EITHER 1 OR 0. 1= plot that parameter | 0= don't plot that
%parameter
%This function is basically the same as DoCoolPlotCombo,  but it allows us
%to control what is being plotted based on conditions. Essentially, it's
%like a toggle based on conditions. 
%BY: Tina BSRI 2021
%[tumor cells, vasc, CART, macrophages, hypoxic cells, antigens?]
%conditions= [1,1,1,0,1,0]; %this array will determine what gets plotted: tumor, vasc, CART,
%macrophages, hypoxia and antigens
conditions=[tumor, vasc, cart, macro, hyp, ant];

% IF het=1, do gradated. if het=0, do binary

iter = iteration;
AA = aa;
AB = ab;
BB = bb;
dirc = strcat('DataGradatedBinary');
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

try
    XYZh=dlmread(strcat('XYZh_t_',int2str(iter),'.txt'));
catch
    XYZh = [-1 -1 -1];
    disp('no hypoxic')
end

%XYZm=dlmread(strcat('XYZm_CART',num2str(AA),'_Trial',num2str(BB),'_t',int2str(iter),'.txt'));
XYZ = dlmread(strcat('XYZloc_t_',int2str(iter),'.txt'));
XYZhreal=XYZh*10;
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
%view([60 30]);
view(250,15);
camlight
lighting gouraud

CoolPlotRadius = 1; % allows smaller than 1 so can see inside...
Coolcellradius = 5;
[spherex,spherey,spherez] = sphere;
spherex = spherex.*CoolPlotRadius.*Coolcellradius;
spherey = spherey.*CoolPlotRadius.*Coolcellradius;
spherez = spherez.*CoolPlotRadius.*Coolcellradius;
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

if conditions(1)==1
    %If we have NO HYPOXIA
    if conditions(5)==0
        [m,n]=size(stemXYZ);
        randsStem = (-5-5).*rand(m,n) + 5;
        stemXYZ=stemXYZ+randsStem;
        [k,t]=size(progXYZ);
        randsProg = (-5-5).*rand(k,t) + 5;
        progXYZ=progXYZ+randsProg;
        disp(strcat("plotting ",num2str(length(XYZ))," cells"));
        %plot stems
        scatter3( stemXYZ(:,1),stemXYZ(:,2),stemXYZ(:,3),40,[0.6350, 0.0780, 0.1840] ,'filled','MarkerEdgeColor','k');
        %plot progenitors
        scatter3( progXYZ(:,1),progXYZ(:,2),progXYZ(:,3),40,[0, 0.75, 0.75] ,'filled','MarkerEdgeColor','k');
        
    else
        %IF we are plotting BOTH hypoxia AND cells:
        %pseudocode:
        %nothyp=HypCell==0 // not hypoxic
        nothypindex=HypCell==0;
        nothyp=XYZreal(nothypindex,:); %CHANGE THIS AFTER RESOLVING THE HYPOXIA ISSUE
        %plot onlyhyp purple
        scatter3(XYZhreal(:,1),XYZhreal(:,2),XYZhreal(:,3),40,[0.4940, 0.1840, 0.5560],'filled','MarkerEdgeColor','k')
        %plot nothyp blue (?)
        scatter3( nothyp(:,1),nothyp(:,2),nothyp(:,3),40,[0, 0.75, 0.75],'filled','MarkerEdgeColor','k');
        %MIXTURE
        
    end
end


XYZkt2 = XYZkt(:,:);
XYZkt2 = XYZkt2*10;
%PLOT CART CELLS: Make them 1/4 of cancer cells
if conditions(3)==1
    %add a random value (-5,5) to each coordinate to see depth on the plot.
    [a,b]=size(XYZkt);
    r = (-5-5).*rand(a,b) + 5;
    disp(strcat("plotting ",num2str(length(XYZkt))," kt cells"));
    scatter3(XYZkt2(:,1)+r(:,1),XYZkt2(:,2)+r(:,2),XYZkt2(:,3)+r(:,3),20, 'green ','filled','MarkerEdgeColor','k')
end

if conditions(4)==1
    %PLOT MACROPHAGES
    scatter3(XYZm(:,1),XYZm(:,2),XYZm(:,3),40, [0.75, 0.75, 0],'filled','MarkerEdgeColor','k')
    
end

if conditions(5)==1
    %Plot hypoxia only if we are not plotting tumor cells. This is because
    %We already plot the case where both tumor and hyp =1. (first for loop)
    if(conditions(1)==0)
scatter3(XYZhreal(:,1),XYZhreal(:,2),XYZhreal(:,3),40,[0.4940, 0.1840, 0.5560],'filled','MarkerEdgeColor','k')
    end
    %PLOT HYPOXIA
end

if conditions(6)==1
    
    [a1,b1]=size(XYZreal);
    r1 = (-5-5).*rand(a1,b1) + 5;
    AntigenHetPlot = dlmread(strcat('AntigenExp_t_',int2str(iter),'.txt'));
    AnitgenHetPlot = round(AntigenHetPlot,3);
    %AntigentHetPlot = [zeros(1,length(AntigenHetPlot));AntigenHetPlot;zeros(1,length(AntigenHetPlot))]; %gradate along green
    colormap(autumn(length(AntigenHetPlot)))
    scatter3(XYZreal(:,1)+r1(:,1),XYZreal(:,2)+r1(:,2),XYZreal(:,3)+r1(:,3),40,AnitgenHetPlot,'filled','MarkerEdgeColor','k');
    colorbar
end
cd ..
%CHANGE THE NAME ACCORDINGLY!!!
%filename = ['PlotEv_AntCartPS_', 'AA',num2str(aa), '_AB', num2str(ab),'_BB', num2str(bb),'_t',num2str(iter), '.tiff'];
filename='shielding_newAngle';
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
%view(150,36);

cd ../../../..
fighandle=f;
end