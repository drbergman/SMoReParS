clearvars;
% Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
% Data point for 5 hours taken out, since it is incommensurate.

tt = [0   10   24   36   48   72  ]';    % hours
tt = tt/24;                       % days

%% Control data
count(:,1) = [0.899 1.340 1.633 2.408 3.557 5.583]';  % thousands of cells
count_std(:,1) = [0.099 0.193 0.207 0.298 0.168 0.364]';  % thousands of cells

%% normalize counts and error to start at 100
temp = 100./count(1,:);
count = count.*temp;
count_std = count_std.*temp;

clear temp

save("data/ExperimentalData.mat",'-v7.3')
