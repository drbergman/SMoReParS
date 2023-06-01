ExperimentalData.mat:
    MAT file that contains experimental data in proper form used for analysis.
    Created by DataFile.m.

OptimalParameters.mat:
    Optimal parameters for SM on cohort_2305311216.
    lambda, alpha, K, a, b are all fixed.
    
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

