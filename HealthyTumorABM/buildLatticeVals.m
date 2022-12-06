function val = buildLatticeVals()

%% Set up occupancy array
val.healthy = 1; % vals to be put in L for healthy cells
val.tumor = 2; % tumor cell val to be put in L for tumor cells
val.healthy_apop = 3; % val to be put in L for apoptotic healthy cells
val.tumor_apop = 4; % val to be put in L for apoptotic tumor cells
val.unused = 5; % val to be put in L for free spaces in computeQ for indexing purposes
