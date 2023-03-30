%% Dose response curve data

Surv  = [100  79.83 65.49 52.59 24.54 10.44  5.57  7.15  5.50]';  % Surviving cell % after 72 hr
Savg = mean(Surv);
OxPt = [0.0001 0.10  0.29  1.00  3.00  9.75 28.65 96.14 300.00]';  % Oxalipltin dose in uM
Surv = [100  79.83 65.49 52.59 24.54 10.44  5.57  7.15  5.50]'/Savg;  % Surviving cell % after 72 hr
Serr = [0    3.98  2.65  4.01  3.53  2.68  4.24  1.28  3.67]'/Savg;




%% 0.75 uM Oxaliplatin
% Cell Cycle Distribution
C1S075 = [90.382 80.918 82.403 85.015 86.974 88.252]';       % % of cells in G0/G1/S

%% 7.55 uM Oxaliplatin
% Cell Cycle Distribution
C1S755 = [90.054 82.258 82.750 80.800 78.313 78.546]';       % % of cells in G0/G1/S

