% this is the original data file Harsh sent.

% data renormalization
% Control data
Cell  = [0.899 1.340 1.633 2.408 3.557 5.583]';  % thousands of cells
Cell075 = [0.899 1.077 1.658 2.059 2.584 3.387]';  % thousands of cells
Cell755 = [0.899 0.960 1.290 1.226 1.254 1.126]';  % thousands of cells
C2M075 = [ 9.629 19.092 17.648 15.006 13.051 11.774]'; % % of cells in G2/M
C2M755 = [ 9.957 17.766 17.256 19.218 21.689 21.475]'; % % of cells in G2/M
Surv  = [100  79.83 65.49 52.59 24.54 10.44  5.57  7.15  5.50]';  % Surviving cell % after 72 hr

C075 = mean(Cell075);
C755 = mean(Cell755);
M075 = mean(C2M075);
M755 = mean(C2M755);
Savg = mean(Surv);

% Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
% Data point for 5 hours taken out, since it is incommensurate.

Time = [0   10   24   36   48   72  ]';    % hours
Time = Time/24;                       % days

% Control data
Cell = [0.899 1.340 1.633 2.408 3.557 5.583]';  % thousands of cells
Cerr = [0.099 0.193 0.207 0.298 0.168 0.364]';  % thousands of cells

% 0.75 uM Oxaliplatin
Cell075 = [0.899 1.077 1.658 2.059 2.584 3.387]'/C075;  % thousands of cells
Cerr075 = [0.104 0.113 0.105 0.370 0.570 0.323]'/C075;  % thousands of cells

% 7.55 uM Oxaliplatin
Cell755 = [0.899 0.960 1.290 1.226 1.254 1.126]'/C755;  % thousands of cells
Cerr755 = [0.098 0.081 0.151 0.248 0.043 0.120]'/C755;  % thousands of cells

% 0.75 uM Oxaliplatin Cell Cycle Distribution
C1S075 = [90.382 80.918 82.403 85.015 86.974 88.252]';       % % of cells in G0/G1/S
C2M075 = [ 9.629 19.092 17.648 15.006 13.051 11.774]'/M075;    % % of cells in G2/M
CerrM075 = [2.3500  1.2090  1.3950  1.2420  1.1640  1.6030]'; %

% 7.55 uM Oxaliplatin Cell Cycle Distribution
C1S755 = [90.054 82.258 82.750 80.800 78.313 78.546]';       % % of cells in G0/G1/S
C2M755 = [ 9.957 17.766 17.256 19.218 21.689 21.475]'/M755;    % % of cells in G2/M
CerrM755 = [1.2550  1.8070  1.8240  2.3920  2.9690  4.4850]';

% Dose response curve data
OxPt = [0.0001 0.10  0.29  1.00  3.00  9.75 28.65 96.14 300.00]';  % Oxalipltin dose in uM
Surv = [100  79.83 65.49 52.59 24.54 10.44  5.57  7.15  5.50]'/Savg;  % Surviving cell % after 72 hr
Serr = [0    3.98  2.65  4.01  3.53  2.68  4.24  1.28  3.67]'/Savg;

% Reshape data
tdata = [Time;Time;Time;Time;OxPt];
cdata = [Cell075;Cell755;C2M075;C2M755;Surv];

data = [tdata cdata];

sigmasq = (mean([Cerr075;Cerr755;(CerrM075/M075);(CerrM755/M755);Serr]))^2
threshold = sigmasq*chi2inv(0.95,5)

