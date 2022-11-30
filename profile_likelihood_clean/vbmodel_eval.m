% Van der ... model

function y = vbmodel_eval(pset, x0, times)
            
    options = odeset('AbsTol',1e-6,'RelTol',1e-8);
    [~,y] = ode15s(@(t,x) vbmodel(x,pset),times,x0,options); % solve ODE systems
        
end
