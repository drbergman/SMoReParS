clearvars;

v = [.01;50];
b = 4;
C = [.75;7.55];

xb = diff(v)/(v(1)/C(1)^b-v(2)/C(2)^b);
blogx = log(xb);
logx = blogx/b;
x = exp(logx);

delta0 = v(1)*((x/C(1))^b+1);

figureOnRight;
scatter(C,v)
hold on;
cc = linspace(0,10,1001);
plot(cc,delta0./((x./cc).^b+1))

%%
bmin = diff(log(v))/diff(log(C));
