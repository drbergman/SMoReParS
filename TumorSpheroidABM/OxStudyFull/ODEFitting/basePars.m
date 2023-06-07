function [p,lb,ub] = basePars(model_type,fixed_pars)


switch model_type
    case "LogisticModel"
        %         1     2   3   4      5      6    7   8       9  10  11
        % p = [lambda;alpha;K;alphaR;alphaP;kalpha;a;delta0;kdelta;b;rho0];

        lambda = 1.534094115133087; lambda_L = 0; lambda_R = 10;
        alpha = 3.299253454440576; alpha_L = 0; alpha_R = 20;
        K = 1029.641812196762; K_L = 400; K_R = 10000;
        alphaR = 1; alphaR_L = 0; alphaR_R = 10;
        alphaP = 1; alphaP_L = 0; alphaP_R = 20;
        kalpha = 1; kalpha_L = 0; kalpha_R = 10;
        a = 4; a_L = 0; a_R = 10;
        delta0 = 1; delta0_L = 0; delta0_R = 50;
        kdelta = 10; kdelta_L = 0; kdelta_R = 30;
        b = 8; b_L = 0; b_R = 10;
        rho0 = .06; rho0_L = 0; rho0_R = 10;

        p = [lambda;alpha;K;alphaR;alphaP;kalpha;a;delta0;kdelta;b;rho0];
        lb = [lambda_L;alpha_L;K_L;alphaR_L;alphaP_L;kalpha_L;a_L;delta0_L;kdelta_L;b_L;rho0_L];
        ub = [lambda_R;alpha_R;K_R;alphaR_R;alphaP_R;kalpha_R;a_R;delta0_R;kdelta_R;b_R;rho0_R];

    case "LogisticModelSimplified"
        %         1     2   3   4      5      6    7        8            9         10
        % p = [lambda;alpha;K;alphaR;alphaP;kalpha;a;low_dose_apop;delta_dose_apop;rho0];

        lambda = 1.534094115133087; lambda_L = 0; lambda_R = 10;
        alpha = 3.299253454440576; alpha_L = 0; alpha_R = 20;
        K = 1029.641812196762; K_L = 400; K_R = 10000;
        alphaR = 1; alphaR_L = 0; alphaR_R = 10;
        alphaP = 1; alphaP_L = 0; alphaP_R = 20;
        kalpha = 1; kalpha_L = 0; kalpha_R = 10;
        a = 4; a_L = 0; a_R = 10;
        low_dose_apop = 1; low_dose_apop_L = 0; low_dose_apop_R = 50;
        delta_dose_apop = 1; delta_dose_apop_L = 0; delta_dose_apop_R = 50;
        rho0 = 0; rho0_L = 0; rho0_R = 10; % in the simplified version, the base value for rho0 is 0 (so if fixed, it is fixed to 0)

        p = [lambda;alpha;K;alphaR;alphaP;kalpha;a;low_dose_apop;delta_dose_apop;rho0];
        lb = [lambda_L;alpha_L;K_L;alphaR_L;alphaP_L;kalpha_L;a_L;low_dose_apop_L;delta_dose_apop_L;rho0_L];
        ub = [lambda_R;alpha_R;K_R;alphaR_R;alphaP_R;kalpha_R;a_R;low_dose_apop_R;delta_dose_apop_R;rho0_R];

    case "HillModel"
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

    otherwise
        error("%s is not a model type that has been specified.",model_type)

end
if nargin>1
    D = parameterOrdering(model_type);
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
