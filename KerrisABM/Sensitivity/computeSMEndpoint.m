function out = computeSMEndpoint(p,~,~)

p(2) = 1 - 1/p(2); % p(2) is nu coming in and theta going forward
opts = odeset("NonNegative",1,"Events",@smEvents);
sol = ode45(@(t,x) odefn(x,p),[0,75],100,opts);

if isempty(sol.ie)
    out = log(sol.y(end));
    return;
end

%% otherwise deal with the events
out = (1/(1-p(2))) * (log(p(1))-log(p(3))); % set it to the (log of) the equilibrium value
