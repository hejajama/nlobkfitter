
#ifndef _DDIS_H
#define _DDIS_H

#include <tuple>

std::tuple<double, double, double, double, double, double> coodinate_inner_products(double x01, double x02, double conj_x01, double conj_x02, double phix0102, double conj_phix0102, double theta20);

double I_ddis_lo_qqbar_L(double Q, double beta, double z, double x01, double conj_x01) ;
double I_ddis_lo_qqbar_T(double Q, double beta, double z, double x01, double conj_x01) ;
double I_ddis_nlo_qqbarg_T_largeM(double Q, double z, double x01sq) ;
double I_ddis_nlo_qqbarg_T_largeQsq(double Q, double beta, double z, double ksq, double r, double rbar) ;
double I_ddis_nlo_qqbarg_L(double Q, double beta, double t, double z0, double z1, double z2, double x01sq, double x02sq, double x21sq, double conj_x01sq, double conj_x02sq, double conj_x21sq) ;
double I_ddis_nlo_qqbarg_L_D3(double Q, double beta, double z0, double z1, double z2, double x01sq, double x02sq, double x21sq, double conj_x01sq, double conj_x02sq, double conj_x21sq) ;
double I_ddis_nlo_qqbarg_T(double Q, double beta, double t, double z0, double z1, double z2, double x01sq, double x02sq, double x21sq, double conj_x01sq, double conj_x02sq, double conj_x21sq) ;
double I_ddis_nlo_qqbarg_T_D3(double Q, double beta, double z0, double z1, double z2, double x01sq, double x02sq, double x21sq, double conj_x01sq, double conj_x02sq, double conj_x21sq) ;

#endif
