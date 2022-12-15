clearvars;

syms w r l

A = [-l-1/w,2*r;1/w,-r-l];

[eV,ev] = eig(A);

V = eV(:,2);
p_lim = V(1)/sum(V);