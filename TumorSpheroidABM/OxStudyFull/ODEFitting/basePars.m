function [p,lb,ub] = basePars(fixed_pars)

%         1     2   3   4      5      6    7   8       9  10  11
% p = [lambda;alpha;K;alphaR;alphaP;kalpha;a;delta0;kdelta;b;rho0];

lambda = 1.534094115133087; lambda_L = 0; lambda_R = 10;
alpha = 3.299253454440576; alpha_L = 0; alpha_R = 20;
K = 1029.641812196762; K_L = 400; K_R = 10000;
alphaR = 1; alphaR_L = 0; alphaR_R = 10;
alphaP = 1; alphaP_L = 0; alphaP_R = 10;
kalpha = 1; kalpha_L = 0; kalpha_R = 10;
a = 4; a_L = 0; a_R = 10;
delta0 = 1; delta0_L = 0; delta0_R = 1000;
kdelta = 10; kdelta_L = 0; kdelta_R = 30;
b = 8; b_L = 0; b_R = 10;
rho0 = .06; rho0_L = 0; rho0_R = 10;

p = [lambda;alpha;K;alphaR;alphaP;kalpha;a;delta0;kdelta;b;rho0];
lb = [lambda_L;alpha_L;K_L;alphaR_L;alphaP_L;kalpha_L;a_L;delta0_L;kdelta_L;b_L;rho0_L];
ub = [lambda_R;alpha_R;K_R;alphaR_R;alphaP_R;kalpha_R;a_R;delta0_R;kdelta_R;b_R;rho0_R];

if nargin>0
    D = parameterOrdering("LogisticModel");
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
