clearvars;
% Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
% Data point for 5 hours taken out, since it is incommensurate.

tt = [0   10   24   36   48   72  ]';    % hours
tt = tt/24;                       % days

doses = [0;0.75;7.55]; % doses in uM

%% Control data
count(:,1) = [0.899 1.340 1.633 2.408 3.557 5.583]';  % thousands of cells
count_std(:,1) = [0.099 0.193 0.207 0.298 0.168 0.364]';  % thousands of cells

% Cell Cycle Distribution (none given)
state2_prop(:,1) = NaN(size(tt));
state2_prop_std(:,1) = NaN(size(tt));

%% 0.75 uM Oxaliplatin
count(:,2) = [0.899 1.077 1.658 2.059 2.584 3.387]';  % thousands of cells
count_std(:,2) = [0.104 0.113 0.105 0.370 0.570 0.323]';  % thousands of cells

% Cell Cycle Distribution
state2_prop(:,2) = [ 9.629 19.092 17.648 15.006 13.051 11.774]' / 100; % proportion of cells in G2/M
state2_prop_std(:,2) = [2.3500  1.2090  1.3950  1.2420  1.1640  1.6030]' / 100; % standard error in proportion of cells in G2/M

%% 7.55 uM Oxaliplatin
count(:,3) = [0.899 0.960 1.290 1.226 1.254 1.126]';  % thousands of cells
count_std(:,3) = [0.098 0.081 0.151 0.248 0.043 0.120]';  % thousands of cells

% Cell Cycle Distribution
state2_prop(:,3) = [ 9.957 17.766 17.256 19.218 21.689 21.475]' / 100; % proportion of cells in G2/M
state2_prop_std(:,3) = [1.2550  1.8070  1.8240  2.3920  2.9690  4.4850]' / 100; % standard error in proportion of cells in G2/M

%% normalize counts and error to start at 100
temp = 100./count(1,:);
count = count.*temp;
count_std = count_std.*temp;

clear temp

save("data/ExperimentalData.mat",'-v7.3')
