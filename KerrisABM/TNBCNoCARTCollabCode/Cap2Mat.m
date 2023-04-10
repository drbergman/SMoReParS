AgentRe = zeros(500,500,500);

nonempty = find(~cellfun('isempty',CapMatrix));
for ii = 1:length(nonempty)
    disp(ii)
    
    A = CapMatrix{ii,1}.ActiveList;
    S = CapMatrix{ii,1}.SegmentList;
    nonempty2 = find(~cellfun('isempty',S));
    for jj = 1:length(nonempty2)
        %disp(strcat(num2str(ii)," ",num2str(jj)));
        Seg = S{jj,1};
        active = A(jj,1)|ii<=2;
        if active
            AgentRe = SegmentFill(AgentRe,Seg.Node1,Seg.Node2,2);
            %disp("a");
        else
            AgentRe = SegmentFill(AgentRe,Seg.Node1,Seg.Node2,1);
        end
    end
end