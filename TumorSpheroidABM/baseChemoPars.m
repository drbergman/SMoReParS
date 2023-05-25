function chemo_pars = baseChemoPars()

chemo_pars.concentration = 0; % in uM
chemo_pars.arrest_function = "linear";

% DNA checks
chemo_pars.dna_check_g1 = false;
chemo_pars.dna_check_s = false;
chemo_pars.dna_check_g2 = false;
chemo_pars.dna_check_m = false;
chemo_pars.dna_check_g1a = false;
chemo_pars.dna_check_sa = false;
chemo_pars.dna_check_g2a = false;
chemo_pars.dna_check_ma = false;

% linear factors connecting concentration and apoptosis probability
chemo_pars.arrest_coeff_g1 = 0.01;
chemo_pars.arrest_coeff_s = 0.01;
chemo_pars.arrest_coeff_g2 = 0.01;
chemo_pars.arrest_coeff_m = 0.01;
chemo_pars.arrest_coeff_g1a = 0.01;
chemo_pars.arrest_coeff_sa = 0.01;
chemo_pars.arrest_coeff_g2a = 0.01;
chemo_pars.arrest_coeff_ma = 0.01;

chemo_pars.arrest_g1_ec50 = 3; % if using chemo_pars.arrest_function = "hill", this is the EC50 for that hill function (Hill coefficient is currently fixed at n=1 for simplicity)
chemo_pars.arrest_s_ec50 = 3; % if using chemo_pars.arrest_function = "hill", this is the EC50 for that hill function (Hill coefficient is currently fixed at n=1 for simplicity)
chemo_pars.arrest_g2_ec50 = 3; % if using chemo_pars.arrest_function = "hill", this is the EC50 for that hill function (Hill coefficient is currently fixed at n=1 for simplicity)
chemo_pars.arrest_m_ec50 = 3; % if using chemo_pars.arrest_function = "hill", this is the EC50 for that hill function (Hill coefficient is currently fixed at n=1 for simplicity)
chemo_pars.arrest_g1a_ec50 = 3; % if using chemo_pars.arrest_function = "hill", this is the EC50 for that hill function (Hill coefficient is currently fixed at n=1 for simplicity)
chemo_pars.arrest_sa_ec50 = 3; % if using chemo_pars.arrest_function = "hill", this is the EC50 for that hill function (Hill coefficient is currently fixed at n=1 for simplicity)
chemo_pars.arrest_g2a_ec50 = 3; % if using chemo_pars.arrest_function = "hill", this is the EC50 for that hill function (Hill coefficient is currently fixed at n=1 for simplicity)
chemo_pars.arrest_ma_ec50 = 3; % if using chemo_pars.arrest_function = "hill", this is the EC50 for that hill function (Hill coefficient is currently fixed at n=1 for simplicity)

% chemo-induced apoptosis-rates
chemo_pars.apoptosis_function = "none";
chemo_pars.apop_c0 = 0; % constant term for chemo-induced increased to apoptosis rate
chemo_pars.apop_c1 = 0; % linear term for chemo-induced increased to apoptosis rate (actual or in Taylor series)
chemo_pars.apop_ec50 = Inf; % if using hill-type function, the ec50 thereof