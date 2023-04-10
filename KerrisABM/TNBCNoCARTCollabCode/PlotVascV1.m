%For plotting vasculature
if mod(time, 120) == 0
%        
%         %mkdir('Run_Images'); V5 now in Control
%         Current_Address = pwd;
% %         figure
% %             hold on
% %             for i = 1: length(CapillaryList)
% %                 PlotNetwork(CapillaryList(i))
% %             end
% % 
% %             SaveAddress = strcat(Current_Address,'/Run_Images/Sticksh', num2str(c(4)), 'm', num2str(c(5)), 't', num2str(time));
% %             %h', num2str(c(3)), 'm', num2str(c(4)),
% %             saveas(gcf, SaveAddress, 'fig'); 
% %             hold off
% 
        figure
            hold on
             p=patch(isosurface(voxelgrid.Agent==1,0));
             set(p,'facecolor','red' ,'edgecolor', 'none');
                daspect([1 1 1])
              % isonormals(voxelgrid.Agent==1,p)
                view([60 30]);

                camlight

                lighting gouraud 
           % SaveAddress = strcat(Current_Address,'/Run_Images/3Dh', num2str(c(4)), 'm', num2str(c(5)), 't', num2str(time));
           % saveas(gcf, SaveAddress, 'fig');    
         hold off
       %  close all
        SaveFig = strcat('Figp', num2str(prolifval), 'm', num2str(migmulti), 'r', num2str(runnum), 't', num2str(time), '.tif');
       print ('-dtiff', '-r300', SaveFig)
%          
%         %Save the txt files 
%         cd Run_images
%          SaveAddress = strcat('Agenth', num2str(c(4)), 'm', num2str(c(5)), 's', num2str(c(6), 2), 't', num2str(time));
%          AgentMatrix = voxelgrid.Agent;
%             save(SaveAddress, 'AgentMatrix');    
%             cd ..
    end