function D = makeMOATDistributions(par_names)

D = dictionary();
for i = 1:numel(par_names)
    switch par_names(i)
        case "carrying_capacity"
            D("carrying_capacity") = makedist("Uniform",500,1e3);
            dmin = 0; dmax = Inf;

        case "g1_to_s"
            D("g1_to_s") = makedist("Normal",24/11,24/110);
            dmin = 1e-5; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        case "s_to_g2"
            D("s_to_g2") = makedist("Normal",3,0.3);
            dmin = 1e-5; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        case "g2_to_m"
            D("g2_to_m") = makedist("Normal",6,0.6);
            dmin = 1e-5; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        case "m_to_g1"
            D("m_to_g1") = makedist("Normal",24,2.4);
            dmin = 1e-5; dmax = Inf; % pick a small, positive minimum transition rate so the initialization can happen even if the icdf is 0

        case "arrest_prob_g1"
            D("arrest_prob_g1") = makedist("Normal",0.05,0.1);
            dmin = 0; dmax = 1;

        case "arrest_prob_g2"
            D("arrest_prob_g2") = makedist("Normal",0.05,0.1);
            dmin = 0; dmax = 1;

        case "apop_rate"
            D("apop_rate") = makedist("Normal",0.1,0.2);
            dmin = 0; dmax = Inf;

        case "move_rate_microns"
            D("move_rate_microns") = makedist("Normal",20,20);
            dmin = 0; dmax = Inf;

        case "occmax_2d"
            D("occmax_2d") = makedist("Uniform",4,8); % cannot make a discrete uniform distribution object, so the draw from this will be floored to produce on integer
            dmin = 0; dmax = 8; % when drawing from this distribution, the algorithm will enforce that occmax_2d<8

        otherwise
            error("have not defined a distribution for %s.",par_names(i))

    end

    D(par_names(i)) = truncate(D(par_names(i)),dmin,dmax);


end
