function p = basePars()

% this loads up what I will consider the base parameter values for the SM
% This SM is given by x' = alpha * x^theta - beta * x

p = zeros(3,1);

p(1) = 1; % alpha
p(2) = 2; % nu (theta = 1 - 1/nu) (this allows nu to be on the open interval (0,inf) which hopefully helps with optimization)
p(3) = 1; % beta