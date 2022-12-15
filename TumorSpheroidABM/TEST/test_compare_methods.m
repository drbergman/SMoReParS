clearvars;

rb = [0.1,3];
Tb = [1/24,16/24];

nr = 12;
nT = 15;

r = linspace(rb(1),rb(2),nr);
T = linspace(Tb(1),Tb(2),nT);

syms I

for ri = 1:nr
    for ti = 1:nT
        a = -1/T(ti);
        b = r(ri)-1/T(ti);
        c = 2*r(ri);

        p = [a,b,c];
        rh = roots(p);
        rh = rh(rh>=0);

        P(ri,ti) = rh/(1+rh);

        S = solve(2==(I+1)*exp(r(ri)*I*T(ti)),I);
        Q(ri,ti) = double(S);


    end
end

%%
figure;
hold on;
mesh(r,T,P');
mesh(r,T,1-Q');
xlabel('r')
ylabel('T')

%%
figure;
hold on;
mesh(r,T,P' - (1-Q'));
xlabel('r')
ylabel('T')