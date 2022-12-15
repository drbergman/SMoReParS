clearvars;

r = 2;
T = 9/24;

syms I

S = solve(2==(I+1)*exp(r*I*T),I);

w1 = linspace(0,T,1001);
w2 = linspace(T,10,10001);

D = r*S*(S+1)*exp(r*(S+1)*T);
C = D*exp(-r*T);

p1 = C*exp(-r*S*w1);
p2 = D*exp(-r*(S+1)*w2);

figure;
plot([w1,w2],[p1,p2])