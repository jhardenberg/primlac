! Common biological parameters

! Equation for light
       double k_bg, k_S, k_D,mN1_D, k_P1, k_P2,I_top

! Equation for Detritus
       double mN1_P1,mN1_P2,l_P1,l_P2
       double   mN1_Z,l_Z,e_z,time_D,mu_Z,P_c

! Equation for Nutrient
       double  mu_P1,N_P1_c,I_P1_c,mu_P2,N_P2_c

! Equation for Phyto 1 and Phyto 2
       double alpha,I_P2_c

! Equation for Zooplankton
       double  gamma_c,Z_c,Fish

       common /light/k_bg, k_S, k_D,mN1_D, k_P1, k_P2,mN1_P1,mN1_P2,l_P1,l_P2,mN1_Z,l_Z,e_z,time_D,mu_z,P_c,I_top
       common /NutrientBiology/ mu_P1,N_P1_c,I_P1_c,mu_P2,N_P2_c
       common /PhytoBiology/alpha,I_P2_c
       common /ZooBiology/gamma_c,Z_c,Fish


