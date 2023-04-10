c = 0;
for x = 1:50
    for y = 1: 50
        for z = 1:50
            if Agentmat(x,y,z)==4
                c = c+1;
            end
        end
    end
end
