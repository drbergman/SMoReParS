function chemo_pars = baseChemoPars()

chemo_pars.concentration = 0; % in uM
chemo_pars.arrest_function = "linear";

% DNA checks
chemo_pars.dna_check_g1 = false;
chemo_pars.dna_check_s = false;
chemo_pars.dna_check_g2 = false;
chemo_pars.dna_check_m = false;
chemo_pars.dna_check_arrest = false;

% linear factors connecting concentration and apoptosis probability
chemo_pars.arrest_coeff_g1 = 0.01;
chemo_pars.arrest_coeff_s = 0.01;
chemo_pars.arrest_coeff_g2 = 0.01;
chemo_pars.arrest_coeff_m = 0.01;
chemo_pars.arrest_coeff_arrest = 0.01;

chemo_pars.arrest_ec50 = 3; % if using chemo_pars.arrest_function = "hill", this is the EC50 for that hill function (Hill coefficient is currently fixed at n=1 for simplicity)
