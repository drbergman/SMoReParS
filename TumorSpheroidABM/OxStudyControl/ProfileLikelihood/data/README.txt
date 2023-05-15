ABMParamEstimates.mat:
    ABM parameters admitted by SMoRe ParS.
    Used just the best fit pars for the ODE.
    See SMoReParS_ProfileSample.m.

ABMParamEstimates_FromProfile.mat:
    ABM parameters admitted by SMoRe ParS.
    Used just the profile from lambda and alpha and a 5^7 grid.
    See SMoReParS_ProfileSample.m.

ABMParamEstimates_FromProfile2.mat:
    ABM parameters admitted by SMoRe ParS.
    Used the profile from lambda and alpha and forced K>597.
    See SMoReParS_ProfileSample.m.

ABMParamEstimates_FromProfile3.mat:
    ABM parameters admitted by SMoRe ParS.
    Used the profile from lambda and alpha on a 3^7 grid.
    See SMoReParS_ProfileSample.m.

ABMParamEstimates_FromProfile_WithK.mat:
    ABM parameters admitted by SMoRe ParS.
    Used full combination information, including on K, to admit ABM parameters.
    SMoReParS_ProfileSample_WithK.m.

LambdaAlphaFn.mat:
    Function approximation for the combination that lambda and alpha are in.
    alpha = f(lambda).
    See FitLambdaAlphaCurve.m for creation of this function.

MultiDimProfileLikelihoods.mat:
    Multidimensional profile of SM from ABM data.
    See MultiDimProfileSMFromABM.m.

ProfileLikelihoods.mat:
    Profile of the SM from the data.
    Profiles are cleaned.
    See ProfileSMFromABM.m.

ProfileLikelihoods_Data.mat:
    Profile of the SM from the data using suboptimal bounds.
    Use ProfileLikelihoods_DataRestricted.mat instead.
    See ProfileSMFromData.m.

ProfileLikelihoods_DataRestricted.mat:
    Profile of the SM from the data using better bounds than ProfileLikelihoods_Data.mat.
    See ProfileSMFromData.m.

