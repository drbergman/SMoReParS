clearvars;

load("ProfileLikelihoods_DataRestricted.mat","out")

xx = out{1}(1,:);
yy = out{1}(2,:);

ind = yy>9.8;
xx(ind) = [];
yy(ind) = [];
f = @(x,p) p(4) + p(1)./((x-p(2)).^p(3));
fobj = @(p) sum((p(1)./((xx-p(2)).^p(3))+p(4)-yy).^2,'all');
p = [1;0;1;1]; % initial guess for y = p(1)/((x-p(2))^p(3)) + p(4);

p = fmincon(fobj,p,[],[],[],[],[0;-Inf;0;-Inf],[Inf;2;Inf;Inf]);

figure;
hold on;
plot(xx,yy)
plot(xx,f(xx,p),"--","LineWidth",2)