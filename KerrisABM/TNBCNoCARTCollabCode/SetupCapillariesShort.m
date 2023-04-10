%SetupCapillariesShortEnd
%This setups the initial capillaries
%This uses fillAgentGrid2
%Last User Kerri: June 22, 2012
%Version 1

%Setup Capillary
C1 = capillary;
C2 = capillary;

%Input Initial Capillaries
            capillary1 = {[30,30,40],[35,27,40];[35,27,40],[35,25,80];[35,25,80],[32,25,120];[32,25,120],[36,30,160];[36,30,160],[32,28,200]};
            capillary2 = {[65,40,40],[65,42,40];[65,42,40],[63,42,80];[63,42,80],[62,41,120];[62,41,120],[64,40,160];[64,40,160],[65,40,200]};
            for i=1:length(capillary1)
                                newSeg = segmentS;
                                newSeg.Node1 = capillary1{i,1};
                                newSeg.Node2 = capillary1{i,2};
                                newSeg.Radius = 2;
                                newSeg.Listpos = i;
                                newSeg.CapillaryNum = 1;
                                newSeg.updateDirection();
                                %seg_list{i} = newSeg;
                                
                                C1.SegmentList{i} = newSeg; 
                                
                                %% Add new cell location to Agent List %ADDED BY KERRI V1.2
                                %Find new cells nodes and radius
                                 Node1pt = newSeg.Node1;
                                 Node2pt = newSeg.Node2;
                                 segradius = newSeg.Radius;
                                %Use FillAgentGrid %NEED TO FIX THIS
                                [voxelgrid] = fillAgentGrid2(voxelgrid, newSeg, pixelsize);

                                %SegmentList{i} = newSeg;
                                %obj.Agents.AssignSpecies(newSeg.Node1(1),newSeg.Node1(2),newSeg.Node1(3),'Segment',newSeg);
                                %obj.Agents.AssignSpecies(newSeg.Node2(1),newSeg.Node2(2),newSeg.Node2(3),'Segment',newSeg);
            end
            
             for j=1:length(capillary2)
                                newSeg = segmentS;
                                newSeg.Node1 = capillary2{j,1};
                                newSeg.Node2 = capillary2{j,2};
                                newSeg.Radius = 2;
                                newSeg.Listpos = j;
                                newSeg.CapillaryNum = 2;
                                newSeg.updateDirection();
                                %seg_list{i} = newSeg;
                                
                                %Setup Capillary
                                %C2 = capillary;
                                C2.SegmentList{j} = newSeg; 
                                
                                %% Add new cell location to Agent List %ADDED BY KERRI V1.2
                                %Find new cells nodes and radius
                                 Node1pt = newSeg.Node1;
                                 Node2pt = newSeg.Node2;
                                 segradius = newSeg.Radius;
                                %Use FillAgentGrid %NEED TO FIX THIS
                                [voxelgrid] = fillAgentGrid2(voxelgrid, newSeg, pixelsize);

                                %SegmentList{i} = newSeg;
                                %obj.Agents.AssignSpecies(newSeg.Node1(1),newSeg.Node1(2),newSeg.Node1(3),'Segment',newSeg);
                                %obj.Agents.AssignSpecies(newSeg.Node2(1),newSeg.Node2(2),newSeg.Node2(3),'Segment',newSeg);
            end


%             %Create two capillary
% C1 = capillary;
% C1.SegmentList{1} = S1; 
% %S1 is now activated as well
% C1.ActiveList = S1.activate(C1.ActiveList);
% 
% C2 = capillary;
% C2.SegmentList{1} = S2; 
% %S1 is now activated as well
% C2.ActiveList = S2.activate(C2.ActiveList);
%             
%             j=1;
%             for i=length(capillary1)+1:length(capillary1)+length(capillary2)
%                                 newSeg = segmentS;
%                                 newSeg.Node1 = capillary2{j,1};
%                                 newSeg.Node2 = capillary2{j,2};
%                                 newSeg.Listpos = i;
%                                 %newSeg.CapillaryNum = 2;
%                                 %seg_list{i} = newSeg;
%                                 SegmentList{i} = newSeg;
%                                 j=j+1;
%                                 %obj.Agents.AssignSpecies(newSeg.Node1(1),newSeg.Node1(2),newSeg.Node1(3),'Segment',newSeg);
%                                 %obj.Agents.AssignSpecies(newSeg.Node2(1),newSeg.Node2(2),newSeg.Node2(3),'Segment',newSeg);
%             end