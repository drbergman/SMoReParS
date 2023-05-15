ExperimentalData.mat:
    Experimental data from old workflow.
    Created by DataFile.m

ODEFittoData.mat:
    Optimal SM parameters fit to data.
    Scales data to start at 100 cells.
    Created by FitODEtoData.

OptimalParameters.mat:
    Optimal SM parameters fit to ABM ran with apoptosis nonzero and varied (cohort_230123120914506 / cohort_230124110240678).
    SM included phase-dependent death terms (I'm pretty sure).
    We are not using this anymore.

OptimalParameters_noapop.mat:
    Optimal SM parameters fit to ABM ran with 0 apoptosis.
    From cohort_230124175743017.
    SM had no apoptosis term.
    This is the main one created at ICERM.

OptimalParameters_Using_DependentStatesAndTimeSeries:
    Optimal SM parameters when assuming dependencies between the state variables and the time series.
    From cohort_230124175743017.
    From same cohort as OptimalParameters_noapop.

OptimalParameters_Using_optimizeSMParsFromABM:
    Optimal SM parameters fit to ABM ran with 0 apoptosis.
    From cohort_230124175743017.
    SM had no apoptosis term.
    This should be the main one moving forward.
    
    