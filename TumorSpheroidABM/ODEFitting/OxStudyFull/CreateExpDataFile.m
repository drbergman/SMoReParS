clearvars;

tt = [0      10     24     36     48     72    ]';        % hours
tt = tt/24;                                             % days

% Control data
data{1} =     [0.899  1.340  1.633  2.408  3.557  5.583]';   % millions of cells
data_std{1} = [0.099  0.193  0.207  0.298  0.168  0.364]';   % millions of cells
% dose = 0.75 uM
prop_in_phase1 = [0.91;0.81;0.82;0.85;0.88;0.9];
data{2} =    [0.899    1     1.633  2.310  2.558  3.289]' .* (prop_in_phase1.*[1,-1] + [0,1]);   % millions of cells
data_std{2} = [0.099  0.193  0.207  0.440  0.850  0.440]' .* (prop_in_phase1.*[1,-1] + [0,1]);   % millions of cells
% dose = 7.55 uM
prop_in_phase1 = [0.91;0.81;0.82;0.80;0.79;0.79];
data{3} =    [0.899   0.920  1.150  1.000  1.000  0.960]' .* (prop_in_phase1.*[1,-1] + [0,1]);   % millions of cells
data_std{3} = [0.099  0.263  0.207  0.210  0.200  0.300]' .* (prop_in_phase1.*[1,-1] + [0,1]);   % millions of cells

%% scale data
factor = zeros(3,1);
for i = 1:3
    factor(i) = 100./sum(data{i}(1,:));
    data_std{i} = data_std{i} .* factor(i);
    data{i} = data{i} .* factor(i);
end

doses = [0;0.75;7.55];

save("data/ExperimentalData.mat","tt","data","data_std","doses","factor","-v7.3")