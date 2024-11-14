% A script to test drawing from ABM parameter space

addpath("~/Documents/MATLAB/myfunctions/") % replace with path (rel or abs) to myfunctions
clearvars

b1 = 2;
n = 101;
u = linspace(0,1,n);
n = length(u);
du = u(2)-u(1);

%% n = 1
A1 = -log(1-(1-exp(-b1))*u)./b1;
if exp(-b1)<eps()
    A1(end) = 1;
end
figureOnRight;
plot(A1,u,"Marker",'.')
axis([0 1 0 1])
% ax_lim = axis;
% if ax_lim(2)>1
%     xlim([0 1])
% end
% if ax_lim(4)>1
%     ylim([0 1])
% end

%% n = 2 and b2 = 0
LHS = @(x) b1*(1-exp(-b1*x)) + exp(-b1*x)*(b1*x+1)-1;
RHS = @(x,v) v*(b1+exp(-b1)-1);
A1 = zeros(n,1);
for i = 1:n
    A1(i) = fsolve(@(x) LHS(x)-RHS(x,u(i)),0.5);
end
figureOnRight;
plot(A1,u,"Marker",'.')
axis([0 1 0 1])

%% n = 2
A1 = zeros(n,1);
b2 = 1e-6;
g = b2-b1;
LHS = @(x) -g*exp(-b1*x) - b1*exp(g*x-b2);
J = @(x) b1*g*(exp(-b1*x)-exp(g*x-b2));
RHS = u*(g*(1-exp(-b1)) + b1*(exp(-b2)-exp(-b1))) - g - b1*exp(-b2);
K = g*(1-exp(-b1)) + b1*(exp(-b2)-exp(-b1));
f = @(x,R) disperse([LHS(x)-R,J(x)]);
% f = @(x,u) (b1-b2)*(exp(-b1*x)-1) + b1*(exp(-b2) - exp((b2-b1)*x-b2)) - u*((1-exp(-b1))*(b2-b1)+b1*(exp(-b2)-exp(-b1)));
% f = @(x,u) (b1-b2)*(exp(b1*x)*(1-exp(b1*(1-x))*u)-1+u) + b1*(exp(-b2) - exp((b2-b1)*x-b2))  - u*b1*(exp(-b2)-exp(-b1));
opts = optimoptions("fsolve","Display","iter-detailed","MaxFunctionEvaluations",1e4,"FunctionTolerance",1e-20,"OptimalityTolerance",1e-16,'SpecifyObjectiveGradient',true,"StepTolerance",1e-9);%,'OutputFcn', @outfun);
% figureOnRight;
for i = 1:n
    % clf;
    if i==1
        A0 = 0;
    else
        A0 = min(1,max(0,A1(i-1,1)-K*du));
    end
    R = RHS(i);
    A1(i) = fsolve(@(x) f(x,R),A0,opts);
end

A2 = -log(1-(1-exp(-b2*A1)).*u)/b2;

%%
figureOnRight;
plot(A1,u,"Marker",".")
axis([0 1 0 1])

%%
hold on;
for i = 1:n
    plot3(A1(i)*ones(n,1),u(i)*ones(n,1),A2(i,:))
end
xlabel('A_1')
ylabel('CDF')
zlabel('A_2')
view(3)

%%
figureOnRight; hold on;
for i = 1:n
    plot(A2(i,:),u)
end

function stop = outfun(x, optimValues, state)
stop = false;
hold on;
plot(x,optimValues.fval,'.');
drawnow
end
