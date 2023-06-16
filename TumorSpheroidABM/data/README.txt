cohort_221217222613078-cohort_230124110240678:
    Cohorts run up to and during first day of ICERM.
    Not using these moving forward.
    These should be suited for deletion.

cohort_230124175743017:
    This is the main cohort for the OxControlStudy.
    The cohort that we made at ICERM and we first performed the full global sensitivity pipeline on.
    This is based on the control case of the oxaliplatin study.

cohort_2303231625:
    This is based on all of the chemotherapeutic concentrations of the oxaliplatin study, including control.
    Big values of arrest_coeff_g1/g2.
    Linked the arrest values instead of varying independently.

cohort_2303271034:
    This is based on all of the chemotherapeutic concentrations of the oxaliplatin study, including control.
    This differs from cohort_2303231625 in that the chemo coefficients are all nonzero: [0.025;0.05;0.075;0.1]

cohort_2303271138:
    This is based on all of the chemotherapeutic concentrations of the oxaliplatin study, including control.
    This differs from cohort_2303271034 in that it also grabbed the coefficient = 0 case from cohort_2303231625.
    That is, cohort_2303271138 = cohort_2303231625 U cohort_2303271034 (U = set union).
    That makes deletion of cohort_2303231625 and cohort_2303271034 ok since they can quickly be recovered from cohort_2303271138 if that subset of chemo coefficients is desired.

cohort_2303301105:
    This is the main cohort we are using for OxStudyFull (the full chemo model).
    This varied the two arrest coefficients independently over a subset of the values used in cohort_2303271138, specifically [0.025;0.05;0.075].
    Thus, it does not fully subsume all the sims in cohort_2303271138.

cohort_2304292140:
    1000 simulations with the same parameter values to test the independence of G1/S and G2/M compartments in the ABM.
    Turns out they were pretty independent (diagonal covariance matrix) by later time points, which I was surpised by.

cohort_2305270925:
    This is a new proposed main cohort for the OxStudyFull (the full chemo model).
    This added arrest compartments for each phase (though only some can become arrested when attempting to transition).
    The arrest function is now a Hill function with rate coefficient being phase-specific, EC50 being phase-independent (+1), and Hill coefficient being fixed at 1.
    Apoptosis occurs only in the arrested compartments in a phase-independent manner.
    The rate follows a Hill function which is 0 at no drug and maxes out at M.chemo_pars.apop_c1 (+2).
    Recovery is phase-independent and constant (+1).
    The above additions introduce 4 new parameters.
    To keep this computationally feasible, we removed the 4 least sensitive parameters from the control study: the 4 transition rates along the cell cycle.

cohort_2305311216:
    This is similar to cohort_2305270925.
    It fixes all ABM parameters varied in the control model, not just transition rates.
    That is, it also fixes carrying capacity, contact inhibition, and migration rate.
    These values are fixed by the average of admitted ABM parameters from the control model (this was not the case for the transition rates in cohort_2305270925).

cohort_2306062212:
    This is similar to cohort_2305311216.
    It fixes the recovery rate to 0 for all arrested compartments.

cohort_2306160948:
    This is similar to cohort_2306062212.
    It sampled one of the accepted ABM parameter vectors from the control study to set the control ABM parameters.
    It selected I=2049==>[3,2,3,1,2,3,3] (carrying capacity, contact inhibition, migration rate, g1->s, s->g2, g2->m, m->g1)
    It fixes the recovery rate to 0 for all arrested compartments.
