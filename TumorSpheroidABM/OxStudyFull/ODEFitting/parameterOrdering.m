function D = parameterOrdering(model_type)

switch model_type
    case "HillModel"
        D = containers.Map(["lambda";"alphaRP";"theta";"VT";"V0";"alphaP";"kalpha";"a";"rho0";"delta0";"kdelta";"b";"alphaR"],1:13);
    case "LogisticModel"
        D = containers.Map(["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"delta0";"kdelta";"b";"rho0"],1:11);
    case "LogisticModelSimplified"
        D = containers.Map(["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"low_dose_apop";"delta_dose_apop";"rho0"],1:10);
    otherwise
        error("%s is not a model type that has been specified.",model_type)
end


