function chemo_pars = baseChemoPars()

chemo_pars.concentration = 0; % in uM

% DNA checks
chemo_pars.dna_check_g1 = false;
chemo_pars.dna_check_s = false;
chemo_pars.dna_check_g2 = false;
chemo_pars.dna_check_m = false;

% linear factors connecting concentration and apoptosis probability
chemo_pars.arrest_coeff_g1 = 0.01;
chemo_pars.arrest_coeff_s = 0.01;
chemo_pars.arrest_coeff_g2 = 0.01;
chemo_pars.arrest_coeff_m = 0.01;
