#include <gsl/gsl_sys.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include "nlodis_config.hpp"
#include "nlodissigmar_massiveq.hpp"
#include "nlodissigmar.hpp"


double I_ddis_lo_qqbar_L(double Q, double beta, double z, double x01, double conj_x01) {
    double res;
    double fac = std::pow(z*(1-z),3.0)*Sq(Q);
    double z1z_root = std::sqrt(z*(1-z));
    double Mx = Q * std::sqrt((1-beta)/beta);
    double j0inner = z1z_root * Mx;
    double k0inner = z1z_root * Q;
    double Phi_0 = x01 * gsl_sf_bessel_J0(j0inner*x01) * gsl_sf_bessel_K0(k0inner*x01);
    double conj_Phi_0 = conj_x01 * gsl_sf_bessel_J0(j0inner*conj_x01) * gsl_sf_bessel_K0(k0inner*conj_x01);
    res = Phi_0 * conj_Phi_0;
    return res;
}

double I_ddis_lo_qqbar_T(double Q, double beta, double z, double x01, double conj_x01) {
    double res;
    double z1z = z * (1-z);
    double fac = Sq(z1z * Q) * (Sq(z) + Sq(1-z));
    double z1z_root = std::sqrt(z1z);
    double Mx = Q * std::sqrt((1-beta)/beta);
    double j1inner = z1z_root * Mx;
    double k1inner = z1z_root * Q;
    double Phi_1 = x01 * gsl_sf_bessel_J1(j1inner*x01) * gsl_sf_bessel_K1(k1inner*x01);
    double conj_Phi_1 = conj_x01 * gsl_sf_bessel_J1(j1inner*conj_x01) * gsl_sf_bessel_K1(k1inner*conj_x01);
    res = Phi_1 * conj_Phi_1;
    return res;
}

double I_ddis_nlo_qqbarg_T_largeM(double Q, double z, double x01sq) {
    double res;
    
    return res;
}

double I_ddis_nlo_qqbarg_T_largeQsq(double Q, double beta, double z, double ksq, double r, double rbar) {
    double res;
    
    return res;
}

double I_ddis_nlo_qqbarg_L(double Q, double beta, double t, double z0, double z1, double z2, double x01sq, double x02sq, double x21sq, double conj_x01sq, double conj_x02sq, double conj_x21sq) {
    double res;
    
    return res;
}

double I_ddis_nlo_qqbarg_L_D3(double Q, double beta, double z0, double z1, double z2, double x01sq, double x02sq, double x21sq, double conj_x01sq, double conj_x02sq, double conj_x21sq) {
    double res;
    
    return res;
}

double I_ddis_nlo_qqbarg_T(double Q, double beta, double t, double z0, double z1, double z2, double x01sq, double x02sq, double x21sq, double conj_x01sq, double conj_x02sq, double conj_x21sq) {
    double res;
    
    return res;
}

double I_ddis_nlo_qqbarg_T_D3(double Q, double beta, double z0, double z1, double z2, double x01sq, double x02sq, double x21sq, double conj_x01sq, double conj_x02sq, double conj_x21sq) {
    double res;
    
    return res;
}

