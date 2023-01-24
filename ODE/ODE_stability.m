clearvars;

syms lambda alpha21 K delta1s delta2m

syms N1s N2m

p = 1-(N1s+N2m)/K;
F = [-lambda*N1s + alpha21*(1+p)*N2m - delta1s*N1s;...
     lambda*N1s - alpha21*N2m - delta2m*N2m];

J = [diff(F,N1s),diff(F,N2m)];

J00 = subs(J,N1s,0);
J00 = subs(J00,N2m,0);

