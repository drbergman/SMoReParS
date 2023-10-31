function D = parameterDisplayNameDictionary(model_type)

switch model_type
    case "HillModel"
        par_names =     ["lambda";"alphaRP";"theta";"VT";"V0";"alphaP";"kalpha";"a";"rho0";"delta0";"kdelta";"b";"alphaR"];
        display_names = ["\lambda";"\alpha_{RP}";"\theta";"V_T";"V_0";"\alpha_P";"k_\alpha";"a";"\rho_0";"\delta_0";"k_\delta";"b";"\alpha_R"];
    case "LogisticModel"
        par_names =     ["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"delta0";"kdelta";"b";"rho0"];
        display_names = ["\lambda";"\alpha";"K";"\alpha_R";"\alpha_P";"k_\alpha";"a";"\delta_0";"k_\delta";"b";"\rho_0"];
    case "LogisticModelSimplified"
        par_names =     ["lambda";"alpha";"K";"alphaR";"alphaP";"kalpha";"a";"low_dose_apop";"delta_dose_apop";"rho0";"delta";"kdelta"];
        display_names = ["\lambda";"\alpha";"K";"\beta_{1S}";"\beta_{2M}";"\kappa_\beta";"m";"d";"\Deltad";"\rho_0";"\delta";"\kappa_\delta"];
    otherwise
        error("%s is not a model type that has been specified.",model_type)
end
D = containers.Map(par_names,display_names);
