%Plot3DFig

%This is used to pull in the capillary data and plot the 3D grid.
%you must load in your data first

tortbin = 1; %1 means use tortuosity, 0 means don't

%Setup Grid
plotgrid.Agent = false(500,500,500);
pixelsize = 1;

%Run through all Capillaries and Plot them
caps = CapMatrix(~cellfun('isempty',CapMatrix));

     for ii=1:length(caps)
         display(ii)
         thesegs = CapMatrix{ii}.SegmentList(~cellfun('isempty',CapMatrix{ii}.SegmentList));
         for jj = 1:length(thesegs)
             if tortbin == 1
                plotgrid = fillAgentGridTortuosity(plotgrid, thesegs{jj}, pixelsize);
             else
                 plotgrid = fillAgentGrid2(plotgrid, thesegs{jj}, pixelsize);
             end
         end
     end
     
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