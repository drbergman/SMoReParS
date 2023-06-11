Profiles_SMFromABM.mat:
    Profile of SM from ABM fixing 5 parameters:
        control parameters (lambda, alpha, K)
        Hill coefficients (a,b)
    Not cleaned.

Profiles_SMFromABM_clean.mat:
    Cleaned version of Profiles_SMFromABM.mat.

Profiles_SMFromABM_Fit_b.mat:
    Profile of SM from ABM fixing 6 parameters:
        control parameters (lambda, alpha, K)
        Hill coefficient controlling arrest rate (a)
    Not cleaned.

Profiles_SMFromABM_Fit_b_clean.mat:
    Cleaned version of Profiles_SMFromABM_Fit_b.mat.

Profiles_SMFromABM_FitAll.mat (not created yet!):
    Profile of SM from ABM fitting all SM parameters.
    Not cleaned.

Profiles_SMFromABM_FitAll_clean.mat (not created yet!):
    Cleaned version of Profiles_SMFromABM_FitAll.mat.

Profiles_SMFromData_LMS.mat:
    Profile of SM using LogisticModelSimplified from data fixing rho0=0
    Not cleaned.

Profiles_SMFromData_LMS_clean.mat:
    Cleaned version of Profiles_SMFromData_LMS.mat.

Profiles_SMFromData.mat:
    Profile of SM from data fixing 5 parameters:
        control parameters (lambda, alpha, K)
        Hill coefficients (a,b)
    Not cleaned.

Profiles_SMFromData_clean.mat:
    Cleaned version of Profiles_SMFromData.mat.

Profiles_SMFromData_Fit_b.mat:
    Profile of SM from data fixing 6 parameters:
        control parameters (lambda, alpha, K)
        Hill coefficient controlling arrest rate (a)
    Not cleaned.

Profiles_SMFromData_Fit_b_clean.mat:
    Cleaned version of Profiles_SMFromData_Fit_b.mat.

Graveyard:
    These were all done with the first version for this model.
    Namely, the one where we used the control model with a death term added; no arrested compartments.
    These remain here simply because I am too scared to delete them right now.

    Profiles_SMFromABM_old.mat:
        This used old workflow and appears to have sometimes resulted in poor profiles.
        See out{3,1} for a profile of the carrying capacity that somehow was between 0 and 2.
        Edit: I seem to have been mistaken about out{3,1}.
            I no longer know if there are poor profiles here, but I will nonetheless continue with the profile from the new workflow.
        Not cleaned.
    
    Profiles_SMFromABM.mat:
        The profile of the SM from the ABM for the new workflow.
        Not cleaned.
    
    Profiles_SMFromABM_extended.mat:
        Same as Profiles_SMFromABM.mat but for each profiled parameter that got close to the boundary, the profile was extended to include this boundary.
        Not cleaned.
    
    Profiles_SMFromABM_clean.mat:
        The profile of the SM from the ABM for the new workflow.
        Cleaned.
        This is exactly Profiles_SMFromABM_extended.mat but cleaned.
        Keeping the original (Profiles_SMFromABM.mat) just in case.
    
    Profiles_SMFromData.mat:
        The profile of the SM from the data for the new workflow.
        Not cleaned.
    
    Profiles_SMFromData_extended.mat:
        Same as Profiles_SMFromData.mat but for each profiled parameter that got close to the boundary, the profile was extended to include this boundary.
        Not cleaned.
    
    Profiles_SMFromData_clean.mat:
        The profile of the SM from the data for the new workflow.
        Cleaned.
        This is exactly Profiles_SMFromData_extended.mat but cleaned.
        Keeping the original (Profiles_SMFromData.mat) just in case.
