%Plot3DFig

%This is used to pull in the capillary data and plot the 3D grid.
%you must load in your data first

tortbin = 1; %1 means use tortuosity, 0 means don't

AA = 10;
AB = 0;
BB = [1:6];

for b = BB
    str = strcat('Data/AA_',num2str(AA),'__BB_',num2str(b),'/');
    %Setup Grid
    plotgrid.Agent = false(500,500,500);
    pixelsize = 2;
    
    %Run through all Capillaries and Plot them
    load(strcat(str,'CapListp_KT',num2str(AA),'_Trial1_t300'));
    caps = CapMatrix(~cellfun('isempty',CapMatrix));
    
    for ii=1:length(caps)
        display(ii)
        thesegs = CapMatrix{ii}.SegmentList(~cellfun('isempty',CapMatrix{ii}.SegmentList));
        for jj = 1:length(thesegs)
            if tortbin == 1
                %plotgrid = fillAgentGridTortuosity(plotgrid, thesegs{jj}, pixelsize);
                plotgrid = fillAgentGridTortNodes(plotgrid, thesegs{jj}, pixelsize);
            else
                plotgrid = fillAgentGrid2(plotgrid, thesegs{jj}, pixelsize);
            end
        end
    end
    
    save(strcat(str,"plotgrid.mat"),'plotgrid');
    
    %Plot it
    figure
    hold on
    p=patch(isosurface(plotgrid.Agent==1,0));
    set(p,'facecolor','red' ,'edgecolor', 'none');
    daspect([1 1 1])
    % isonormals(voxelgrid.Agent==1,p)
    view([60 30]);
    
    camlight
    
    lighting gouraud
    % SaveAddress = strcat(Current_Address,'/Run_Images/3Dh', num2str(c(4)), 'm', num2str(c(5)), 't', num2str(time));
    % saveas(gcf, SaveAddress, 'fig');
    hold off
end