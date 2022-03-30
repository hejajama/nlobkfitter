
#ifndef _NLODIS_SIGMAR_MASSIVE_H
#define _NLODIS_SIGMAR_MASSIVE_H


double L_dip( double Q, double z, double mf );
double OmegaL_V( double Q, double z, double mf );
double G_qg_int(int a, int b, double Qbar, double mf, double x2, double x3, double omega, double lambda);
double G_qg(int a, int b, double y_u, double y_t, double Qbar, double mf, double x2, double x3, double omega, double lambda);

double ILdip_massive_LiLogConst(double Q, double z1, double x01sq, double mf) ;
double ILdip_massive_Iab(double Q, double z1, double x01sq, double mf, double xi) ;
double ILdip_massive_Icd(double Q, double z1, double x01sq, double mf, double xi, double x) ;
double ILNLOqg_massive_dipole_part(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
double ILNLOqg_massive_tripole_part(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_u, double y_t) ;
double ILNLOqg_massive_dipole_part_symm(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
double ILNLOqg_massive_tripole_part_symm(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_u, double y_t) ;
double ILNLOqg_massive_dipole_part_unintegrated(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2);
double ILNLOqg_massive_tripole_part_unintegrated(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2);

#endif // _NLODIS_SIGMAR_MASSIVE_H
