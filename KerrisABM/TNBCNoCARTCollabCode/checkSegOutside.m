AgentRe = zeros(500,500,500);

nonempty = find(~cellfun('isempty',CapMatrix));
for ii = 1:length(nonempty)
    disp(ii)
    
    S = CapMatrix{ii,1}.SegmentList;
    nonempty2 = find(~cellfun('isempty',S));
    for jj = 1:length(nonempty2)
        %disp(strcat(num2str(ii)," ",num2str(jj)));
        Seg = S{jj,1};
        Node1 = Seg.Node1;
        Node2 = Seg.Node2;
        N1 = sum(Node1>1000|Node1<1);
        N2 = sum(Node2>1000|Node2<1);
        if(N1>0)
            error(strcat("Error: OOB", Node1));
        end
        if(N2>0)
            error(strcat("Error: OOB", Node2));
        end
    end
end