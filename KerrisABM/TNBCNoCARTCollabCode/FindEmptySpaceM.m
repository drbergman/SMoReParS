function [newposition, CellTypeList] = FindEmptySpaceM(CAgrid, Position, indmat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FindEmptySpace
%This fuction will check the positions of nearby cells and see if there is
%any empty spaces
%CAgrid is the binary grid which has the locations of the cells
%Position is the [x y z] of the current cell
%indmat is the indices of the Matrix
%This records the types of cells next to the cell of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            onfirst = Position == 1;
            onlast = size(CAgrid) == Position; 
            if sum(onfirst >0) || sum(onlast > 0)
                lowx = max(1, Position(1)-1);
                lowy = max(1, Position(2)-1);                
                lowz = max(1, Position(3)-1);
                higx = min(size(CAgrid,1), Position(1)+1);                
                higy = min(size(CAgrid,2), Position(2)+1);                
                higz = min(size(CAgrid,3), Position(3)+1);
                
                PosList = {};
                CellTypeList = {};
                
                for ix = lowx:higx
                    for iy = lowy:higy
                        for iz = lowz:higz
                            thePos = [ix iy iz];
                            if CAgrid(ix,iy,iz) == 0
                                %Add to empty position list
                                PosList = cat(1, PosList, thePos);
                            else
                                CellTypeList = cat(1, CellTypeList, CAgrid(ix,iy,iz));
                            end
                        end
                    end
                end
                
                if ~isempty(PosList)
                    randlist = randperm(length(PosList));
                     rpos = randlist(1);
                     newposition = PosList{rpos,:};
                else
                    %There is no open position
                    newposition = {};
                end
                
            else
                %The cell is not on the edge of the grid
                 Neighborgrid = CAgrid(Position(1)-1:Position(1)+1,Position(2)-1:Position(2)+1, Position(3)-1:Position(3)+1);
                 %newposition = Neighborgrid;
                 %[row,col,len] = ind2sub(size(Neighborgrid),find(Neighborgrid == 0));
                 NoNeigh = find(Neighborgrid == 0);

                 if ~isempty(NoNeigh) %check there is an open position
                     %Choose random direction
                     randlist = randperm(length(NoNeigh));
                     rpos = NoNeigh(randlist(1));
                     thesubs = indmat(rpos,:);

                     rowi = thesubs(1) - 2; %Put row into -1 0 1 
                     coli = thesubs(2) - 2; %Put col into -1 0 1 
                     leni = thesubs(3) - 2; %Put len into -1 0 1 
                     newposition = [Position(1)+rowi, Position(2)+coli, Position(3)+leni];
                 else
                     newposition = {};
                 end
            end % if sum
         end %function
  
    
