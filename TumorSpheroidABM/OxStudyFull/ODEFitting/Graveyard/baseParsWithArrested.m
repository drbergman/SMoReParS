function [p,lb,ub] = baseParsWithArrested(fixed_pars)

lambda = 5.96; lambda_L = 0; lambda_R = 100;
alphaRP = 0.99; alphaRP_L = 0; alphaRP_R = 1000;
theta = 1;  theta_L = 0; theta_R = 100;
VT = 548;  VT_L = 0; VT_R = 10000;
V0 = 55;  V0_L = 0; V0_R = 20000;
alphaP = 8.69;  alphaP_L = 0; alphaP_R = 200;
kalpha = 21.67;  kalpha_L = 0; kalpha_R = 10000;
a = 0.5;  a_L = 0; a_R = 100;
rho0 = 0.06;  rho0_L = 0; rho0_R = 20; % not listed in their Table 1 under Oxaliplatin, using the value reported in Taxol model of Table 1
delta0 = 0.96;  delta0_L = 0; delta0_R = 100;
kdelta = 9.22;  kdelta_L = 0; kdelta_R = 10000;
b = 1;  b_L = 0; b_R = 200;
alphaR = 3.0;  alphaR_L = 0; alphaR_R = 1000;

p = [lambda;alphaRP;theta;VT;V0;alphaP;kalpha;a;rho0;delta0;kdelta;b;alphaR];

lb = [lambda_L;alphaRP_L;theta_L;VT_L;V0_L;alphaP_L;kalpha_L;a_L;rho0_L;delta0_L;kdelta_L;b_L;alphaR_L];
ub = [lambda_R;alphaRP_R;theta_R;VT_R;V0_R;alphaP_R;kalpha_R;a_R;rho0_R;delta0_R;kdelta_R;b_R;alphaR_R];

if nargin>0
    D = parameterOrdering("HillModel");
    for i = 1:numel(fixed_pars)
        if ~isKey(D,fixed_pars(i))
            if fixed_pars(i)==""
                continue
            else
                warning("%s is not a parameter in this model.\n",fixed_pars(i))
            end
        end
        I = D(fixed_pars(i));
        lb(I) = p(I);
        ub(I) = p(I);
    end
end
