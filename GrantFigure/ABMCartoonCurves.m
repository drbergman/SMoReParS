% this script will make the cartoon of ABM growth curves

clearvars;
n = 4; % number of growth curves per parameter point
colors = ["#FBB040";"#2E3192"];

t = [0 1 2 3 4];
[~,ybar(2,:)] = ode45(@(t,x) 2*x*(1-x),t,0.01);
[~,ybar(1,:)] = ode45(@(t,x) 2*x*(1-1.2*x),t,0.01);
ystd = .08./(1+(t-2).^(2));
ystd(1) = 0;
% ybar = [1,2,5,7,8;...
%         1,2.5,5.4,7.5,9];
% ystd = [0,.2,1,.5,.2;...
%         0,0.3,1.1,0.2,.1];

y = ybar + ystd.*randn([2,5,n]);

figure; hold on
for i = 1:2
    plot(t,squeeze(y(i,:,:)),"Color",colors(i),"LineWidth",0.5)
end
xticks([])
yticks([])
ax = gca;
ax.Units = "pixels";
ax.Position(3:4) = 72;
