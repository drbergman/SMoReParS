function p = basePars()

% this loads up what I will consider the base parameter values for the SM
% This SM is given by x' = alpha * x^theta - beta * x

p = zeros(3,1);

p(1) = 1; % alpha
p(2) = 0.5; % theta
p(3) = 1; % beta