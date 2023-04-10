cd Data/AA_5__AB_20/Run_Images
%load('CapListp_t_20.mat');
cd ../../..
nodes = zeros(30,3);
nonempty = find(~cellfun('isempty',CapMatrix));
for ii = 1:length(nonempty)
  nodes(ii, :) = CapMatrix{ii}.SegmentList{1}.Node1;
end
isDup = zeros(length(nodes),1);
for ii = 1:length(nodes)
    isDup(ii) =sum( nodes(ii,1)==nodes(:,1)&nodes(ii,2)==nodes(:,2)&nodes(ii,3)==nodes(:,3));
end

disp([nodes isDup])
isDup = isDup>1;
disp(sum(isDup));