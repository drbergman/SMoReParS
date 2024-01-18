function D = makeSMParameterDistributionsDictionary(par_names)

% the bounds here come from the [10%,90%] quantiles of the parameter
% distributions found in fitting the ODE to the ABM on the lattice
D = dictionary();
for i = 1:numel(par_names)
    switch par_names(i)
        case "lambda"
%             D("lambda") = makedist("Uniform",0.8689,1.5561); % for the fit with apop and all 5 parameters
            D("lambda") = makedist("Uniform",0.9915,1.4684); % for the fit without apop and 3 parameters
            dmin = 0; dmax = Inf;

        case "alpha"
%             D("alpha") = makedist("Uniform",3.5965,6.2461); % for the fit with apop and all 5 parameters
            D("alpha") = makedist("Uniform",3.9991,5.9492); % for the fit without apop and 3 parameters
            dmin = 0; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        case "K"
%             D("K") = makedist("Uniform",488.1592,4032.9); % for the fit with apop and all 5 parameters
            D("K") = makedist("Uniform",526.2866,3143.0); % for the fit without apop and 3 parameters
            dmin = 0; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        case "delta"
            D("delta") = makedist("Uniform",4.0002e-8,0.4017);
            dmin = 0; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        case "g1_prop0"
            D("g1_prop0") = makedist("Uniform",0.7425,0.8322);
            dmin = 0; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        otherwise
            error("have not defined a distribution for %s.",par_names(i))

    end

    D(par_names(i)) = truncate(D(par_names(i)),dmin,dmax);


end
