function cycle_pars = baseCyclePars()

% transition rates from one phase to the next in per day; all based on typical durations of the phases
cycle_pars.g1_to_s = 24/11;
cycle_pars.s_to_g2 = 24/8;
cycle_pars.g2_to_m = 24/4;
cycle_pars.m_to_g1 = 24/1;

% DNA checks
cycle_pars.dna_check_g1 = true;
cycle_pars.dna_check_s = true;
cycle_pars.dna_check_g2 = true;
cycle_pars.dna_check_m = false;

cycle_pars.arrest_prob_g1 = 0.2;
cycle_pars.arrest_prob_s = 0.2;
cycle_pars.arrest_prob_g2 = 0.2;
cycle_pars.arrest_prob_m = 0.2;
