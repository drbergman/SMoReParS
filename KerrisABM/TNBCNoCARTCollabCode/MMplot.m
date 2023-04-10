X = 1:125000;
Y = arrayfun(@helper,X);
divis = arrayfun(@(x) .2,X);
disp(find(divis==Y))
plot(X,Y)
hold on
plot(X,divis)
scatter(25000,helper(25000))
hold off

function prob = helper(numberofcellscart)
    offset = 25000;
    cartRandDeath = 0.051;
    CartBreakpt = 25000;
    halfpoint = 32000;%37500;
    if numberofcellscart < CartBreakpt
     %if false
        prob = cartRandDeath;
    else
        %numberofcellscart = numberofcellscart-offset;
        %MM = (.225-cartRandDeath)*numberofcellscart/(halfpoint-offset + numberofcellscart);
        %MM = (.225-cartRandDeath)*(numberofcellscart-CartBreakpt)/(7000 + numberofcellscart-CartBreakpt);
        MM = (.225-cartRandDeath)*(numberofcellscart-CartBreakpt)/(25000 + numberofcellscart-CartBreakpt);
        %prob = MM;
        prob = MM+0.051;
    end
end