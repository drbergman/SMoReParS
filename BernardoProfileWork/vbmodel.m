% Kamaldeen Okuneye 2021-04-04
function dx = vbmodel(x,p)
   
% K = pset(1); theta = pset(2); mu_P = pset(3);
dx = p(3)*x*((p(1)/x)^(1/p(2)) - 1); 
