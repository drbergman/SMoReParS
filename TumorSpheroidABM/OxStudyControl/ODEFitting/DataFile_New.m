clearvars;
% Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
% Data point for 5 hours taken out, since it is incommensurate.

t = [0   10   24   36   48   72  ]';    % hours
t = t/24;                       % days

cohort_size = 1; % only one cohort represented in this data

%% Control data
D.A = [0.899 1.340 1.633 2.408 3.557 5.583]';  % thousands of cells
D.S = [0.099 0.193 0.207 0.298 0.168 0.364]';  % thousands of cells

%% normalize counts and error to start at 100
temp = 100./D.A(1);
D.A = D.A.*temp;
D.S = D.S.*temp;

C = {[]};

save("data/ExperimentalData_New.mat","t","D","C","cohort_size",'-v7.3')
