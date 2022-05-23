
#ifndef _NLODIS_SIGMAR_MASSIVE_H
#define _NLODIS_SIGMAR_MASSIVE_H


// ------------------- HELPER FUNCTIONS --------------------------

double L_dip( double Q, double z, double mf );
double OmegaL_V( double Q, double z, double mf );
double OmegaT_V_unsymmetric( double Q, double z, double mf );
double OmegaT_N_unsymmetric( double Q, double z, double mf );

double IT_V1_unsymmetric( double Q, double z, double mf, double r, double xi );
double IT_VMS1_unsymmetric( double Q, double z, double mf, double r, double xi );
double IT_V2_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u );
double IT_VMS2_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u );
double IT_N_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u );

double IT_dipole_jk_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_dipole_jkm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);

double IT_tripole_jk_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_tripole_jkm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_tripole_F_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_tripole_Fm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double IT_tripole_jk_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u);
double IT_tripole_jkm_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u);
double IT_tripole_F_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u);
double IT_tripole_Fm_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u);
double IT_tripole_jk_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2);
double IT_tripole_jkm_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2);
double IT_tripole_F_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2);
double IT_tripole_Fm_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2);




double G_qg_int(int a, int b, double Qbar, double mf, double x2, double x3, double omega, double lambda);
double G_qg(int a, int b, double y_u, double y_t, double Qbar, double mf, double x2, double x3, double omega, double lambda);


// ---------------- LONGITUDINAL INTEGRANDS ------------------------

double ILdip_massive_LiLogConst(double Q, double z1, double x01sq, double mf) ;
double ILdip_massive_Iab(double Q, double z1, double x01sq, double mf, double xi) ;
double ILdip_massive_Icd(double Q, double z1, double x01sq, double mf, double xi, double x) ;
double ILNLOqg_massive_dipole_part_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double ILNLOqg_massive_tripole_part_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double ILNLOqg_massive_tripole_part_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u);
double ILNLOqg_massive_tripole_part_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2);

// For debugging purposes:
double ILNLOqg_massive_dipole_part(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
double ILNLOqg_massive_tripole_part(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_u, double y_t) ;
double ILNLOqg_massive_dipole_part_symm(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
double ILNLOqg_massive_tripole_part_symm(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_u, double y_t) ;
double ILNLOqg_massive_dipole_part_unintegrated(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2);
double ILNLOqg_massive_tripole_part_unintegrated(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2);




// ---------------- TRANSVERSE INTEGRANDS ------------------------

double ITdip_massive_0(double Q, double z1, double x01sq, double mf) ;
double ITdip_massive_1(double Q, double z1, double x01sq, double mf, double xi) ;
double ITdip_massive_2(double Q, double z1, double x01sq, double mf, double y_chi, double y_u) ;

double ITNLOqg_massive_dipole_part_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double ITNLOqg_massive_tripole_part_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq);
double ITNLOqg_massive_tripole_part_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u);
double ITNLOqg_massive_tripole_part_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2);


#endif // _NLODIS_SIGMAR_MASSIVE_H
