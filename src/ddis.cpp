#include <math>
#include <vector>
#include <numerics>
#include <algorithm> // for transform
#include <functional> // for plus

#include <gsl/gsl_sys.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include "ddis.hpp"
#include "nlodis_config.hpp"
#include "nlodissigmar.hpp"
#include "nlodissigmar_massiveq.hpp"

std::tuple<double, double, double, double, double, double> coodinate_inner_products(double x01, double x02, double conj_x01, double conj_x02, double phix0102, double conj_phix0102, double theta20){
    std::vector<double> x01v = {x01, 0};
    std::vector<double> x02v = {x02*std::cos(phix0102), x02*std::sin(phix0102)};
    std::vector<double> x21v;
    std::transform(x01v.begin(), x01v.end(), x02v.begin(), std::back_inserter(x21v), std::minus<float>()); // x21v = x01v - x02v
    std::vector<double> conj_x01v = {conj_x01, 0};
    std::vector<double> conj_x02v = {conj_x02*std::cos(conj_phix0102), conj_x02*std::sin(conj_phix0102)};
    std::vector<double> conj_x21v;
    std::transform(conj_x01v.begin(), conj_x01v.end(), conj_x02v.begin(), std::back_inserter(conj_x21v), std::minus<float>());
    double net_rotation = theta20 + phix0102 - conj_phix0102;
    double cosnet = std::cos(net_rotation);
    double sinnet = std::sin(net_rotation);
    conj_x01v = std::vector<double> {conj_x01v[0]*cosnet - conj_x01v[1]*sinnet, conj_x01v[0]*sinnet + conj_x01v[1]*cosnet};
    conj_x02v = std::vector<double> {conj_x02v[0]*cosnet - conj_x02v[1]*sinnet, conj_x02v[0]*sinnet + conj_x02v[1]*cosnet};
    conj_x21v = std::vector<double> {conj_x21v[0]*cosnet - conj_x21v[1]*sinnet, conj_x21v[0]*sinnet + conj_x21v[1]*cosnet};

    double x02_dot_x21 = std::inner_product(x02v.begin(), x02v.end(), x21v.begin(), 0);
    double conj_x02_dot_conj_x21 = std::inner_product(conj_x02v.begin(), conj_x02v.end(), conj_x21v.begin(), 0);
    double conj_x02_dot_x02 = std::inner_product(conj_x02v.begin(), conj_x02v.end(), x02v.begin(), 0);
    double conj_x02_dot_x21 = std::inner_product(conj_x02v.begin(), conj_x02v.end(), x21v.begin(), 0);
    double conj_x21_dot_x02 = std::inner_product(conj_x21v.begin(), conj_x21v.end(), x02v.begin(), 0);
    double conj_x21_dot_x21 = std::inner_product(conj_x21v.begin(), conj_x21v.end(), x21v.begin(), 0);

    std::tuple<double, double, double, double, double, double> dot_products {x02_dot_x21, conj_x02_dot_conj_x21, conj_x02_dot_x02, conj_x02_dot_x21, conj_x21_dot_x02, conj_x21_dot_x21};
    return dot_products;
}

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

