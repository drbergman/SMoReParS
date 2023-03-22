% this script will make the cartoon of SM growth curves

clearvars;
n = 4; % number of growth curves per parameter point
colors = ["#FBB040";"#2E3192"];

t = [0 1 2 3 4];
r = [linspace(1.8,2.2,4);linspace(1.6,1.9,4)];
K = [linspace(0.8,1.1,4);linspace(0.7,1,4)];
K = K(:,randperm(4));
y = zeros(2,length(t),n);
for i = 1:4
    for j = 1:2
        [~,y(j,:,i)] = ode45(@(t,x) r(j,i)*x*(1-x/K(j,i)),t,0.01);
    end
end

y = flip(y,1);

figure; hold on
for i = 1:2
    plot(t,squeeze(y(i,:,:)),"Color",colors(i),"LineWidth",0.5)
end
xticks([])
yticks([])
ax = gca;
ax.Units = "pixels";
ax.Position(3:4) = 72;
