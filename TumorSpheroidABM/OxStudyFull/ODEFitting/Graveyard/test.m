clearvars;

p = basePars();

doses = [0;0.75;7.55];
n = 100;

for i = 1:n
    for j = 1:numel(doses)
        dose = doses(j);

        dose_arrest_factor = 1/(1+(p(6)/dose)^p(7));
        dose_apoptosis_factor = p(8)/(1+(p(9)/dose)^p(10));

        sol = ode23(@(t,x) odefn(x,p,dose_arrest_factor,dose_apoptosis_factor),[0 3],[90;10;0;0]);
    end
end