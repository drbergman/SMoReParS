function dx = odefn(x,p,dose_arrest_factor,dose_apoptosis_factor,recovery_rate)

% the ODE function defining the SM
%         1     2   3   4      5      6    7   8       9  10  11
% p = [lambda,alpha,K,alphaR,alphaP,kalpha,a,delta0,kdelta,b,rho0]

dx = [-p(1) * x(1) + (2 - sum(x)/p(3))*(p(2) * x(2)) - p(4) * dose_arrest_factor * x(1) + recovery_rate * (x(3)+x(4));
     p(1) * x(1) - p(2) * x(2) - p(5) * dose_arrest_factor * x(2);
     p(4) * dose_arrest_factor * x(1) - recovery_rate * x(3) - dose_apoptosis_factor * x(3);
     p(5) * dose_arrest_factor * x(2) - recovery_rate * x(4) - dose_apoptosis_factor * x(4)];

%% Explanation of the ODE
% dx = zeros(4,1); % [G1/S; G2/M; arrested G1/S; arrested G2/M];
% 
% transition_to_proliferating = p(1) * x(1);
% transition_to_resting = p(2) * x(2);
% arrest_of_resting = p(4) * dose_arrest_factor * x(1);
% arrest_of_proliferating = p(5) * dose_arrest_factor * x(2);
% recovery_from_resting_arrested = recovery_rate * x(3);
% recovery_from_proliferating_arrested = recovery_rate * x(4);
% apoptosis_of_resting_arrested = dose_apoptosis_factor * x(3);
% apoptosis_of_proliferating_arrested = dose_apoptosis_factor * x(4);
% 
% dx(1) = -transition_to_proliferating + (2 - sum(x)/p(3))*transition_to_resting - arrest_of_resting + recovery_from_resting_arrested + recovery_from_proliferating_arrested;
% dx(2) = transition_to_proliferating - transition_to_resting - arrest_of_proliferating;
% dx(3) = arrest_of_resting - recovery_from_resting_arrested - apoptosis_of_resting_arrested;
% dx(4) = arrest_of_proliferating - recovery_from_proliferating_arrested - apoptosis_of_proliferating_arrested;
