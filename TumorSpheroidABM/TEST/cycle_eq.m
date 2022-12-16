clearvars;

syms r1 r2 r3 r4

A = sym(zeros(4));
A(1,[1,4]) = [-r1,2*r4];
A(2,[1,2]) = [r1,-r2];
A(3,[2,3]) = [r2,-r3];
A(4,[3,4]) = [r3,-r4];

[eV,ev] = eig(A);