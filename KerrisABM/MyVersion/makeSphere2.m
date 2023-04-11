%makeSphere2
%This fills a discrete grid with sphere of radius R, around the point p1
%This file also stores the indices and segment type for each Agent
%This file changes floor and ceil to round
%Version 2.3
%Last Revision: Kerri May 11, 2012

%% Kerri's version
% function [grid] = makeSphere2(grid, p1, R, indices) 
% gs = size(grid.Agent);
% for i = max(1,round(p1(1)) - R^2):min(gs(1),round(p1(1)) + R^2)
% 
%     for j = max(1,round(p1(2)) - R^2):min(gs(2),round(p1(2)) + R^2)
%         for k = max(1,round(p1(3)) - R^2):min(gs(3),round(p1(3)) + R^2)
%             if ((i - p1(1))^2 + (j - p1(2))^2 + (k-p1(3))^2) <= R^2
%                 grid.Agent(i,j,k) = true;
%                 if int16(indices(1)) == 0
%                     error('CapInd is 0!!!')
%                 end
%                 grid.CapInd(i,j,k) = int16(indices(1));
%                 grid.PosInd(i,j,k) = int16(indices(2));
%                 %tempCell.Type = segtype;
%                 %grid.Cell(i,j,k) = tempCell;
%             end
%         end
%     end
% end

%% Daniel's update for Kerri's version
% function [grid] = makeSphere2(grid, p1, R, indices)
% if int16(indices(1)) == 0
%     error('CapInd is 0!!!')
% end
% gs = size(grid.Agent);
% for i = max(1,round(p1(1)) - R):min(gs(1),round(p1(1)) + R)
%     for j = max(1,round(p1(2)) - R):min(gs(2),round(p1(2)) + R)
%         for k = max(1,round(p1(3)) - R):min(gs(3),round(p1(3)) + R)
%             if ((i - p1(1))^2 + (j - p1(2))^2 + (k-p1(3))^2) <= R^2
%                 grid.Agent(i,j,k) = true;
%                 grid.CapInd(i,j,k) = int16(indices(1));
%                 grid.PosInd(i,j,k) = int16(indices(2));
%                 %tempCell.Type = segtype;
%                 %grid.Cell(i,j,k) = tempCell;
%             end
%         end
%     end
% end

%% Daniel's update for Kerri's version v2
% function [grid] = makeSphere2(grid, p1, R, indices)
% if int16(indices(1)) == 0
%     error('CapInd is 0!!!')
% end
% gs = size(grid.Agent);
% for i = max(1,round(p1(1)) - R):min(gs(1),round(p1(1)) + R)
%     remaining_r_x = R^2 - (i-p1(1))^2;
%     if remaining_r_x < 0
%         continue;
%     end
%     remaining_r_x = sqrt(remaining_r_x);
%     for j = max(1,ceil(p1(2) - remaining_r_x)):min(gs(2),floor(p1(2) + remaining_r_x))
%         remaining_r_y = R^2 - (i-p1(1))^2 - (j-p1(2))^2;
%         if remaining_r_y < 0
%             continue;
%         end
%         remaining_r_y = sqrt(remaining_r_y);
%         for k = max(1,ceil(p1(3) - remaining_r_y)):min(gs(3),floor(p1(3) + remaining_r_y))
%             if ((i - p1(1))^2 + (j - p1(2))^2 + (k-p1(3))^2) <= R^2 % this check really is necessary
%                 grid.Agent(i,j,k) = true;
%                 grid.CapInd(i,j,k) = int16(indices(1));
%                 grid.PosInd(i,j,k) = int16(indices(2));
%                 %tempCell.Type = segtype;
%                 %grid.Cell(i,j,k) = tempCell;
%             end
%         end
%     end
% end

%% Daniel's update for Kerri's version using grid
function [grid] = makeSphere2(grid, p1, R, indices,onlyAgentGridNeeded)
if nargin==4
    onlyAgentGridNeeded = false;
end
if int16(indices(1)) == 0
    error('CapInd is 0!!!')
end
if isfield(grid,"sz")
    gs = grid.sz;
else
    gs = size(grid.Agent);
end
[X,Y,Z] = ndgrid(max(1,ceil(p1(1)-R)):min(gs(1),floor(p1(1)+R)),max(1,ceil(p1(2)-R)):min(gs(2),floor(p1(2)+R)),max(1,ceil(p1(3)-R)):min(gs(3),floor(p1(3)+R)));
D = (X-p1(1)).^2+(Y-p1(2)).^2+(Z-p1(3)).^2;
I = find(D<=R^2);
ind = sub2ind(gs,X(I),Y(I),Z(I));
grid.Agent(ind) = true;
if ~onlyAgentGridNeeded
    grid.CapInd(ind) = int16(indices(1));
    grid.PosInd(ind) = int16(indices(2));
end

%% Daniel's version for integers
% function [grid] = makeSphere2(grid, p1, R, indices)
% if int16(indices(1)) == 0
%     error('CapInd is 0!!!')
% end
% gs = size(grid.Agent);
% for i = 0:R
%     if i==0
%         xi = p1(1);
%         nxi = 1;
%     else
%         xi = p1(1)+i*[-1,1];
%         if p1(1)+i>gs(1)
%             xi(2)=[];
%         end
%         if p1(1)<=i
%             xi(1)=[];
%         end
%         nxi = length(xi);
%     end
%     for j = 0:floor(sqrt(R^2-i^2))
%         if j==0
%             yi = p1(2);
%             nyi = 1;
%         else
%             yi = p1(2)+j*[-1,1];
%             if p1(2)+j>gs(2)
%                 yi(2)=[];
%             end
%             if p1(2)<=j
%                 yi(1)=[];
%             end
%             nyi = length(yi);
%         end
%         for k = 0:floor(sqrt(R^2-i^2-j^2))
%             if k==0
%                 zi = p1(3);
%                 nzi = 1;
%             else
%                 zi = p1(3)+k*[-1,1];
%                 if p1(3)+k>gs(3)
%                     zi(2)=[];
%                 end
%                 if p1(3)<=k
%                     zi(1)=[];
%                 end
%                 nzi = length(zi);
%             end
%             for ii = 1:nxi
%                 for jj = 1:nyi
%                     for kk = 1:nzi
%                         grid.Agent(xi(ii),yi(jj),zi(kk)) = true;
%                         grid.CapInd(xi(ii),yi(jj),zi(kk)) = int16(indices(1));
%                         grid.PosInd(xi(ii),yi(jj),zi(kk)) = int16(indices(2));
%                     end
%                 end
%             end
%         end
%     end
% end