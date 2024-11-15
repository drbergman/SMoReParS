ExperimentalData.mat:
    MAT file that contains experimental data in proper form used for analysis.
    Created by DataFile.m.

SMFitToABM_LMS.mat:
    LMS stands for LogisticModelSimplified.
    Optimal parameters for SM on cohort_2306062212.
    rho0 is fixed at 0.

SMFitToABM_LMS_bounded.mat:
    LMS stands for LogisticModelSimplified.
    Optimal parameters for SM on cohort_2306062212.
    rho0 is fixed at 0.
    The chemo-induced death rate at low dose is bounded above zero.

SMFitToABM_LMS_2.mat:
    LMS stands for LogisticModelSimplified.
    Optimal parameters for SM on cohort_2306160948.
    rho0 is fixed at 0.
    The chemo-induced death rate at low dose is bounded above zero.

SMFitToABM_Fit_b.mat:
    Optimal parameters for SM on cohort_2305311216.
    lambda, alpha, K, a are all fixed.

SMFitToABM_FitAll.mat:
    Similar to SMFitToABM_Fit_b.mat.
    All SM parameters fit.

SMFitToABM.mat:
    Similar to SMFitToABM_Fit_b.mat.
    This also fixes b compared to that one.

SMFitToData_LMS.mat:
    Optimal parameters for SM using the LogisticModelSimplified on experimental data.
    rho0 is fixed at 0. 

SMFitToData_LMS_bounded.mat:
    Optimal parameters for SM using the LogisticModelSimplified on experimental data.
    rho0 is fixed at 0. 
    The chemo-induced death rate at low dose is bounded above zero.

SMFitToData_Fit_b.mat:
    Optimal parameters for SM on experimental data.
    lambda, alpha, K, a are all fixed.

SMFitToData_FitAll.mat:
    Similar to SMFitToData_Fit_b.mat.
    All SM parameters fit.

SMFitToData.mat:
    Similar to SMFitToData_Fit_b.mat.
    This also fixes b compared to that one.
    
Graveyard:
    From the first SM that was the control with chemo-induced death.
    
    ODEFitToData.mat:
        Best fit of SM to data.
        The SM is the control model with arrested compartments added.
    
    OptimalParameters.mat:
        Optimal parameters for SM on cohort_2303271138 assuming that the chemo-induced death rate is constant across the two ODE state variables.
    
    OptimalParameters_phase_dependent_death.mat:
        Optimal parameters for SM on cohort_2303271138 assuming that the chemo-induced death rate depends on the duration spent in each phase.
        Specifically, longer duration in one phase corresponds with slower death rate to agree with the arrest only occuring at the intra-phase transition (G1->S or G2->M).
    
    OptimalParameters_UnLinkedHill:
        Optimal parameters for SM on cohort_2303301105 assuming that the chemo-induced death rate follows a Hill function that is independent for each phase.
        The Hill coefficient was fixed at 3.
        After running the script FitODEtoABM.m, I must have manually deleted this ind from the parameter dimension.

