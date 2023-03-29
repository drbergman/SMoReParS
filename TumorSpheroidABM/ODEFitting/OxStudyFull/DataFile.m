{\rtf1\ansi\ansicpg1252\cocoartf2708
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 ArialMT;\f1\froman\fcharset0 Times-Roman;}
{\colortbl;\red255\green255\blue255;\red71\green71\blue71;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c34902\c34902\c34902;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww28600\viewh16880\viewkind0
\deftab720
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 % data renormalization
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 % Control data
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Cell\'a0 \'a0 = [0.899\'a0 1.340\'a0 1.633\'a0 2.408\'a0 3.557\'a0 5.583]'; \'a0 % thousands of cells
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Cell075 = [0.899\'a0 1.077\'a0 1.658\'a0 2.059\'a0 2.584\'a0 3.387]'; \'a0 % thousands of cells
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Cell755 = [0.899\'a0 0.960\'a0 1.290\'a0 1.226\'a0 1.254\'a0 1.126]'; \'a0 % thousands of cells
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 C2M075\'a0 = [ 9.629\'a0 19.092\'a0 17.648\'a0 15.006\'a0 13.051\'a0 11.774]'; % % of cells in G2/M
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 C2M755\'a0 = [ 9.957\'a0 17.766\'a0 17.256\'a0 19.218\'a0 21.689\'a0 21.475]'; % % of cells in G2/M
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Surv\'a0 \'a0 = [100\'a0 \'a0 79.83\'a0 65.49\'a0 52.59\'a0 24.54\'a0 10.44 \'a0 5.57 \'a0 7.15 \'a0 5.50]'; \'a0 % Surviving cell % after 72 hr
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 C075 = mean(Cell075);
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 C755 = mean(Cell755);
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 M075 = mean(C2M075);
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 M755 = mean(C2M755);
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Savg = mean(Surv);
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 % Data from Jang et al. Cancer Res Treat 2002;34:372. Millions of cells.
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 % Data point for 5 hours taken out, since it is incommensurate.
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 Time = [0\'a0 \'a0 \'a0 10 \'a0 \'a0 24 \'a0 \'a0 36 \'a0 \'a0 48 \'a0 \'a0 72\'a0 \'a0 ]';\'a0 \'a0 \'a0 \'a0 % hours
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Time = Time/24; \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 % days
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 % Control data
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Cell = [0.899\'a0 1.340\'a0 1.633\'a0 2.408\'a0 3.557\'a0 5.583]'; \'a0 % thousands of cells
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Cerr = [0.099\'a0 0.193\'a0 0.207\'a0 0.298\'a0 0.168\'a0 0.364]'; \'a0 % thousands of cells
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 % 0.75 uM Oxaliplatin\'a0
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Cell075 = [0.899\'a0 1.077\'a0 1.658\'a0 2.059\'a0 2.584\'a0 3.387]'/C075; \'a0 % thousands of cells
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Cerr075 = [0.104\'a0 0.113\'a0 0.105\'a0 0.370\'a0 0.570\'a0 0.323]'/C075; \'a0 % thousands of cells
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 % 7.55 uM Oxaliplatin
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Cell755 = [0.899\'a0 0.960\'a0 1.290\'a0 1.226\'a0 1.254\'a0 1.126]'/C755; \'a0 % thousands of cells
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Cerr755 = [0.098\'a0 0.081\'a0 0.151\'a0 0.248\'a0 0.043\'a0 0.120]'/C755; \'a0 % thousands of cells
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 % 0.75 uM Oxaliplatin Cell Cycle Distribution
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 C1S075 = [90.382\'a0 80.918\'a0 82.403\'a0 85.015\'a0 86.974\'a0 88.252]'; \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 % % of cells in G0/G1/S
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 C2M075 = [ 9.629\'a0 19.092\'a0 17.648\'a0 15.006\'a0 13.051\'a0 11.774]'/M075;\'a0 \'a0 \'a0 \'a0 % % of cells in G2/M
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 CerrM075 = [2.3500\'a0 \'a0 1.2090\'a0 \'a0 1.3950\'a0 \'a0 1.2420\'a0 \'a0 1.1640\'a0 \'a0 1.6030]';\'a0 %\'a0
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 % 7.55 uM Oxaliplatin Cell Cycle Distribution
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 C1S755 = [90.054\'a0 82.258\'a0 82.750\'a0 80.800\'a0 78.313\'a0 78.546]'; \'a0 \'a0 \'a0 \'a0 \'a0 \'a0 % % of cells in G0/G1/S
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 C2M755 = [ 9.957\'a0 17.766\'a0 17.256\'a0 19.218\'a0 21.689\'a0 21.475]'/M755;\'a0 \'a0 \'a0 \'a0 % % of cells in G2/M
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 CerrM755 = [1.2550\'a0 \'a0 1.8070\'a0 \'a0 1.8240\'a0 \'a0 2.3920\'a0 \'a0 2.9690\'a0 \'a0 4.4850]';
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 % Dose response curve data
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 OxPt = [0.0001\'a0 0.10 \'a0 0.29 \'a0 1.00 \'a0 3.00 \'a0 9.75\'a0 28.65\'a0 96.14 300.00]'; \'a0 % Oxalipltin dose in uM
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Surv = [100\'a0 \'a0 79.83\'a0 65.49\'a0 52.59\'a0 24.54\'a0 10.44 \'a0 5.57 \'a0 7.15 \'a0 5.50]'/Savg; \'a0 % Surviving cell % after 72 hr
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 Serr = [0 \'a0 \'a0 \'a0 3.98 \'a0 2.65 \'a0 4.01 \'a0 3.53 \'a0 2.68 \'a0 4.24 \'a0 1.28 \'a0 3.67]'/Savg;
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 % Reshape data
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 tdata = [Time;Time;Time;Time;OxPt];
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 cdata = [Cell075;Cell755;C2M075;C2M755;Surv];
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 data = [tdata cdata];
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\pardeftab720\sa320\partightenfactor0

\f0\fs48 \cf2 \strokec2 sigmasq = (mean([Cerr075;Cerr755;(CerrM075/M075);(CerrM755/M755);Serr]))^2
\f1\fs24 \cf0 \strokec3 \

\f0\fs48 \cf2 \strokec2 threshold = sigmasq*chi2inv(0.95,5)
\f1\fs24 \cf0 \strokec3 \
\pard\pardeftab720\partightenfactor0
\cf0 \
}