% this script looked at the combination of (lambda,alpha) that was found in
% main_data.m and fits a curve to that alpha = f(lambda)

clearvars;

load("ProfileLikelihoods_DataRestricted.mat","out")

xx = out{1}(1,:); % lambda
yy = out{1}(2,:); % alpha

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