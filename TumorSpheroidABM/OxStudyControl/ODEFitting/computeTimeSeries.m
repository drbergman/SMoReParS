function out = computeTimeSeries(p,tt,~,fn_opts,data)

% computes the time series solution for the SM at time points tt. Always
% uses initial conditions of [90;10];

if fn_opts.condition_on_previous
    out = zeros(length(tt),2);
    out(1,:) = data(1,:);
    sol = ode45(@(t,x) odefn(x,p),tt(1:2),data(1,:)');
    out(2,:) = deval(sol,tt(2));
    for i = 3:length(tt)
        sol = odextend(sol,[],tt(i),data(i-1,:)');
        out(i,:) = deval(sol,tt(i));
    end
else
    [~,out] = ode45(@(t,x) odefn(x,p),tt,[90;10]);
end
if size(data,2)==1
    out = sum(out,2);
end

% sol = ode45(@(t,x) odefn(x,p),[0 3],[90;10]);
% 
% out = deval(sol,tt)';