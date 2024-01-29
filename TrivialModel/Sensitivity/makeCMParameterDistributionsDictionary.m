function [D,I] = makeCMParameterDistributionsDictionary(par_names, opts)

arguments
    par_names string
    opts.distribution = "Uniform"
end


D = dictionary(); % distribution dictionary
I = dictionary(["a","b"],1:2); % index dictionary
for i = 1:numel(par_names)
    switch par_names(i)
        case "a"
            switch opts.distribution
                case "Uniform"
                    D("a") = makedist("Uniform",1,3);
                case "Normal"
                    D("a") = makedist("Normal",2,0.05);
            end
            dmin = 1; dmax = 3;

        case "b"
            switch opts.distribution
                case "Uniform"
                    D("b") = makedist("Uniform",4,7);
                case "Normal"
                    D("b") = makedist("Normal",5.5,0.05);
            end
            dmin = 4; dmax = 7;

        otherwise
            error("have not defined a distribution for %s.",par_names(i))

    end

    D(par_names(i)) = truncate(D(par_names(i)),dmin,dmax);


end