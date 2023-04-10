%AgentRe = zeros(500,500,500);

nonempty = find(~cellfun('isempty',CapMatrix));
for ii = 1:length(nonempty)
    disp(ii)
    
    S = CapMatrix{ii,1}.SegmentList;
    nonempty2 = find(~cellfun('isempty',S));
    for jj = 1:length(nonempty2)
        %disp(strcat(num2str(ii)," ",num2str(jj)));
        Seg = S{jj,1};
        Node1 = Seg.Node1(1);
        Node2 = Seg.Node1(2);
        Node3 = Seg.Node1(3);
        if (Node1 == 499 && Node2 ==210 && Node3 == 2)
            error(strcat("Cap ",num2str(ii)," Seg ",num2str(jj)));
        end

    end
end