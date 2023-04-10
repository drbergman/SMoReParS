%SetupCapillariesEnd
%This setups the initial capillaries
%This uses fillAgentGrid2
%Last User Kerri: March 12, 2013
%Version 1

%Setup Capillary
C1 = capillary;
C2 = capillary;

%Input Initial Capillaries
            capillary1 = {[30,30,20],[35,27,40];[35,27,40],[35,25,80];[35,25,80],[32,25,120];[32,25,120],[36,30,160];[36,30,160],[32,28,200];[32,28,200],[30,30,240];[30,30,240],[28,30,280];[28,30,280],[27,30,320];[27,30,320],[25,30,360];[25,30,360],[28,30,400]};
            capillary2 = {[465,440,20],[465,442,40];[465,442,40],[463,442,80];[463,442,80],[462,441,120];[462,441,120],[464,440,160];[464,440,160],[465,440,200];[465,440,200],[463,440,240];[463,440,240],[467,442,280];[467,442,280],[466,444,320];[466,444,320],[465,442,360];[465,442,360],[465,440,400]};
            for ik=1:length(capillary1)
                                newSeg = segmentS;
                                newSeg.Node1 = capillary1{ik,1};
                                newSeg.Node2 = capillary1{ik,2};
                                newSeg.Radius = 2;
                                newSeg.Listpos = ik;
                                newSeg.CapillaryNum = 1;
                                newSeg.updateDirection();
                                %seg_list{i} = newSeg;
                                
                                %fill in SegmentMatrix
                                SegmentMatrix(ik,1:3) =  capillary1{ik,1};
                                SegmentMatrix(ik,4:6) =  capillary1{ik,2};
                                SegmentMatrix(ik,7) =  1;
                                SegmentMatrix(ik,8) =  ik;
                                SegmentMatrix(ik,9) =  1;
                                
                                
                                C1.SegmentList{ik} = newSeg; 
                                
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
                                jj = ik+j;
                                %fill in SegmentMatrix
                                SegmentMatrix(jj,1:3) =  capillary1{j,1};
                                SegmentMatrix(jj,4:6) =  capillary1{j,2};
                                SegmentMatrix(jj,7) =  1;
                                SegmentMatrix(jj,8) =  j;
                                SegmentMatrix(jj,9) =  1;
                                
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