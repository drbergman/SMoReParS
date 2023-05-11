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

