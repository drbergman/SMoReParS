function [XYZcart] = GetCartSpawnLoc(matNodes, cartnum,xyz,indmat)
%GetCartSpawnLoc: gets positions for car-t cells to enter the simulation

%   This function takes in the matrix of mature nodes, and the number of
%   new car-t cells that are being created. It returns n random points along
%   the mature capilaries where n = cartnum.
%   Here we assume that a CAR-T cell will emerge from any point in the
%   mature vasculature with equal probability. 

    noZeros = matNodes(all(matNodes,2),:);
    [numpts,~] = size(noZeros);
    idx = randsample(1:numpts , cartnum);
    XYZcart = ceil(noZeros(idx,:)/20);
    XYZcart(ind2sub(size(XYZcart),find(XYZcart==0)))=1; 
    
    XYZcells = num2cell(XYZcart, 2);
    allLocPos = cellfun(@(x) x-2+indmat,XYZcells,'UniformOutput',false);
    XYZcart = cell2mat(allLocPos);
    XYZcart(ind2sub(size(XYZcart),find(XYZcart==0)))=1; 
    XYZcart(ind2sub(size(XYZcart),find(XYZcart>50)))=50; 
    XYZcart = unique(XYZcart,'rows');
    
    XYZcart = setdiff(XYZcart,xyz,'stable','rows');
    
    cartnum = min(cartnum,length(XYZcart));
    finalidx = randsample(1:length(XYZcart),cartnum);
    XYZcart = XYZcart(finalidx,:);
    

end

