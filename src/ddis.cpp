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

// ----- Final state transfer functions
double I2_Mx(double Mx, double z0, double z1, double Y01){
    return 1/(4*M_PI)*gsl_sf_bessel_J0(std::sqrt(z0*z1)*Mx*Y01);
}

double I2_Delta(double t, double coordinates_todo){
    return 1/(4*M_PI)*gsl_sf_bessel_J0(std::sqrt(std::abs(t))*coordinates_todo)
}

double I3_Mx(double Mx, double z0, double z1, double z2, double Y012){
    return 2*(z0*z1*z2)/Sq(4*M_PI)*Mx/Y012*gsl_sf_bessel_J1(Mx*Y012);
}

double I3_Delta(double t, double coordinates_todo){
    return 1/(4*M_PI)*gsl_sf_bessel_J0(std::sqrt(std::abs(t))*coordinates_todo);
}


// ----- LO & APPROX qqbarg DDIS impact factors -----

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
    double bessel_inner = std::sqrt(x01sq*z*(1-z))*Q;
    double psi_LO_sq = 1/M_PI * z*(1-z) * (Sq(2*z - 1) + 1) * Sq(Q*gsl_sf_bessel_K1(bessel_inner));
    double impact_fac_integrand = psi_LO_sq;
    
    return impact_fac_integrand;
}

double I_ddis_nlo_qqbarg_T_largeQsq(double Q, double beta, double z, double ksq, double r, double rbar) {
    double beta_over_z = beta / z;
    double splitting_fun = Sq(beta_over_z) + Sq(1-beta_over_z);
    double j2inner = std::sqrt((1-z)*ksq);
    double k2inner = std::sqrt(z*ksq);
    double impact_fac_integrand = Sq(ksq) * std::log(Sq(Q)/ksq) * splitting_fun
                                  * gsl_sf_bessel_J2(j2inner*r) * gsl_sf_bessel_J2(j2inner*rbar)
                                  * gsl_sf_bessel_K2(k2inner*r) * gsl_sf_bessel_K2(k2inner*rbar); // r * rbar multiplier not here since polar coord jacobians have been kept in nlodissigmar.cpp
    
    return impact_fac_integrand;
}


// 
// ----- NLO DDIS impact factors -----
// 

std::tuple<double, double, double, double, double, double, double> coodinate_inner_products(double x01, double x02, double conj_x01, double conj_x02, double phix0102, double conj_phix0102, double theta20){
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
    double conj_x01_dot_x01 = std::inner_product(conj_x01v.begin(), conj_x01v.end(), x01v.begin(), 0);
    double conj_x02_dot_x02 = std::inner_product(conj_x02v.begin(), conj_x02v.end(), x02v.begin(), 0);
    double conj_x02_dot_x21 = std::inner_product(conj_x02v.begin(), conj_x02v.end(), x21v.begin(), 0);
    double conj_x21_dot_x02 = std::inner_product(conj_x21v.begin(), conj_x21v.end(), x02v.begin(), 0);
    double conj_x21_dot_x21 = std::inner_product(conj_x21v.begin(), conj_x21v.end(), x21v.begin(), 0);

    std::tuple<double, double, double, double, double, double, double> dot_products {x02_dot_x21, conj_x02_dot_conj_x21, conj_x01_dot_x01, conj_x02_dot_x02, conj_x02_dot_x21, conj_x21_dot_x02, conj_x21_dot_x21};
    return dot_products;
}

double I_ddis_nlo_qqbarg_L(double Q, double beta, double t, double z0, double z1, double z2, double x01, double x02, double phix0102, double conj_x01, double conj_x02, double conj_phix0102, double theta20) {
    // compute t dependent factor and fall back to FL_D(3)
    double coordinates_todo = 0;
    double t_dep = I3_Delta(t, coordinates_todo);
    double FLD3 = I_ddis_nlo_qqbarg_L_D3(Q, beta, z0, z1, z2, x01, x02, phix0102, conj_x01, conj_x02, conj_phix0102, theta20);
    return t_dep*FLD3;
}

double I_ddis_nlo_qqbarg_L_D3(double Q, double beta, double z0, double z1, double z2, double x01, double x02, double phix0102, double conj_x01, double conj_x02, double conj_phix0102, double theta20) {
    double res;
    double x01sq = Sq(x01);
    double x02sq = Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double conj_x01sq = Sq(conj_x01);
    double conj_x02sq = Sq(conj_x02);
    double conj_x21sq=conj_x01sq+conj_x02sq-2.0*sqrt(conj_x01sq*conj_x02sq)*cos(conj_phix0102);
    std::tuple<double, double, double, double, double, double> dot_products = coodinate_inner_products(x01, x02, conj_x01, conj_x02, phix0102, conj_phix0102, theta20);
    double conj_x01_dot_x01 = std::get<2>(dot_products); // needed for Y012
    double conj_x02_dot_x02 = std::get<3>(dot_products);
    double conj_x02_dot_x21 = std::get<4>(dot_products);
    double conj_x21_dot_x02 = std::get<5>(dot_products);
    double conj_x21_dot_x21 = std::get<6>(dot_products);

    double Mx = Q*std::sqrt((1-beta)/beta);
    double Y012 = z0*z1*(conj_x01sq + x01sq - 2 * conj_x01_dot_x01)
                + z1*z2*(conj_x21sq + x21sq - 2 * conj_x21_dot_x21)
                + z0*z2*(conj_x02sq + x02sq - 2 * conj_x02_dot_x02);
    double I_Mx = I3_Mx(Mx, z0, z1, z2, Y012);

    double besk_inner_QX = Q*std::sqrt(z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq);
    double conj_besk_inner_QX = Q*std::sqrt(z1*(1.0 - z1 - z2)*conj_x01sq + z2*(1.0 - z1 - z2)*conj_x02sq + z2*z1*conj_x21sq);
    double BesK0 = gsl_sf_bessel_K0(besk_inner_QX);
    double conj_BesK0 = gsl_sf_bessel_K0(conj_besk_inner_QX);
    double big_brace_terms = 2 * (
        Sq(z1) * (
            (2 * z0 * (z0 + z2) + Sq(z2)) * (conj_x02_dot_x02/(conj_x02sq*x02sq) - (-conj_x21_dot_x02)/(2*conj_x21sq*x02sq) - (-conj_x02_dot_x21)/(2*conj_x02sq*x21sq))
            + Sq(z2)/2 * (-conj_x02_dot_x21/(conj_x02sq*x21sq) + (-conj_x21_dot_x02)/(conj_x21sq*x02sq))
            )
    ); // factor of 2 from combining the terms under the q<->qbar exchange symmetry
    double impact_fac_integrand = 4 * z0 * z1 * Sq(Q) * BesK0 * conj_BesK0 * big_brace_terms;

    res = I_Mx * impact_fac_integrand / (z0 * z1 * z2); // Remember to divide by the z_i under the integration measure
    return res;
}

double I_ddis_nlo_qqbarg_T(double Q, double beta, double t, double z0, double z1, double z2, double x01, double x02, double phix0102, double conj_x01, double conj_x02, double conj_phix0102, double theta20) {
    // compute t dependent factor and fall back to FT_D(3)
    double coordinates_todo = 0;
    double t_dep = I3_Delta(t, coordinates_todo);
    double FTD3 = I_ddis_nlo_qqbarg_T_D3(Q, beta, z0, z1, z2, x01, x02, phix0102, conj_x01, conj_x02, conj_phix0102, theta20);
    return t_dep*FTD3;
}

double I_ddis_nlo_qqbarg_T_D3(double Q, double beta, double z0, double z1, double z2, double x01, double x02, double phix0102, double conj_x01, double conj_x02, double conj_phix0102, double theta20) {
    double res;
    double x01sq = Sq(x01);
    double x02sq = Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double conj_x01sq = Sq(conj_x01);
    double conj_x02sq = Sq(conj_x02);
    double conj_x21sq=conj_x01sq+conj_x02sq-2.0*sqrt(conj_x01sq*conj_x02sq)*cos(conj_phix0102);
    std::tuple<double, double, double, double, double, double> dot_products = coodinate_inner_products(x01, x02, conj_x01, conj_x02, phix0102, conj_phix0102, theta20);
    double x02_dot_x21 = std::get<0>(dot_products); // not in L
    double x21_dot_x02 = x02_dot_x21;
    double conj_x02_dot_conj_x21 = std::get<1>(dot_products); // not in L
    double conj_x21_dot_conj_x02 = conj_x02_dot_conj_x21;
    double conj_x01_dot_x01 = std::get<2>(dot_products); // needed for Y012
    double conj_x02_dot_x02 = std::get<3>(dot_products);
    double conj_x02_dot_x21 = std::get<4>(dot_products);
    double conj_x21_dot_x02 = std::get<5>(dot_products);
    double conj_x21_dot_x21 = std::get<6>(dot_products);
    // 3-index x_i+j;k inner products -- sc for semicolon, p for 'plus'
    // b^2 term
    double conj_x_0p2sc1_dot_x_0p2sc1 = Sq(z0/(z0+z2))*conj_x02_dot_x02 - z0/(z0+z2)*(-conj_x02_dot_x21 - conj_x21_dot_x02) + conj_x21_dot_x21;
    double conj_x_0p2sc1_dot_conj_x20 = -z0/(z0+z2)*conj_x02sq + (-conj_x21_dot_conj_x02); // note x20 is used here, and sign flip is accounted for (instead of x02=x0-x2)
    double x_0p2sc1_dot_x20 = -z0/(z0+z2)*x02sq + (-x21_dot_x02); // note x20 is used here, and sign flip is accounted for (instead of x02=x0-x2)
    double conj_x_0p2sc1_dot_x20 = -z0/(z0+z2)*conj_x02_dot_x02 + (-conj_x21_dot_x02);
    double x_0p2sc1_dot_conj_x20 = -z0/(z0+z2)*conj_x02_dot_x02 + (-conj_x02_dot_x21);
    // c^2 term
    double conj_x_0sc1p2_dot_x_0sc1p2 = Sq(z1/(z1+z2))*conj_x21_dot_x21 - z1/(z1+z2)*(-conj_x02_dot_x21 - conj_x21_dot_x02) + conj_x02_dot_x02;
    double conj_x_0sc1p2_dot_conj_x21 = -(-conj_x02_dot_conj_x21) + z1/(z1+z2)*conj_x21sq;
    double x_0sc1p2_dot_x21 = -(-x02_dot_x21) + z1/(z1+z2)*x21sq;
    double conj_x_0sc1p2_dot_x21 = -(-conj_x02_dot_x21) + z1/(z1+z2)*conj_x21_dot_x21;
    double x_0sc1p2_dot_conj_x21 = -(-conj_x21_dot_x02) + z1/(z1+z2)*conj_x21_dot_x21;
    // interference term
    double conj_x_0p2sc1_dot_x_0sc1p2 = z0/(z0+z2)*conj_x02_dot_x02 - z0*z1/((z0+z2)*(z1+z2))*(-conj_x02_dot_x21) - (-conj_x21_dot_x02) + z1/(z1+z2)*conj_x21_dot_x21;
    double conj_x_0sc1p2_dot_x_0p2sc1 = z0/(z0+z2)*conj_x02_dot_x02 - (-conj_x02_dot_x21) - z0*z1/((z0+z2)*(z1+z2))*(-conj_x21_dot_x02) + z1/(z1+z2)*conj_x21_dot_x21;
    // double conj_x_0p2sc1_dot_conj_x20 = ; // already in b^2 terms
    double conj_x_0p2sc1_dot_x21 = -z0/(z0+z2)*(-conj_x02_dot_x21) + conj_x21_dot_x21;
    // double x_0sc1p2_dot_x21 = ; // already in c^2
    double x_0sc1p2_dot_conj_x20 = - conj_x02_dot_x02 + z1/(z1+z2)*(-conj_x02_dot_x21);
    // double conj_x_0sc1p2_dot_conj_x21 = ; already in c^2
    double conj_x_0sc1p2_dot_x20 = - conj_x02_dot_x02 + z1/(z1+z2)*(-conj_x21_dot_x02);
    // double x_0p2sc1_dot_x20 = ; // already in b^2
    double x_0p2sc1_dot_conj_x21 = -z0/(z0+z2)*(-conj_x21_dot_x02) + conj_x21_dot_x21;

    double Mx = Q*std::sqrt((1-beta)/beta);
    double Y012 = z0*z1*(conj_x01sq + x01sq - 2 * conj_x01_dot_x01)
                + z1*z2*(conj_x21sq + x21sq - 2 * conj_x21_dot_x21)
                + z0*z2*(conj_x02sq + x02sq - 2 * conj_x02_dot_x02);
    double I_Mx = I3_Mx(Mx, z0, z1, z2, Y012);

    double X012 = std::sqrt(z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq);
    double conj_X012 = std::sqrt(z1*(1.0 - z1 - z2)*conj_x01sq + z2*(1.0 - z1 - z2)*conj_x02sq + z2*z1*conj_x21sq);
    double besk_inner_QX = Q*X012;
    double conj_besk_inner_QX = Q*conj_X012;
    double BesK1 = gsl_sf_bessel_K1(besk_inner_QX)/X012;
    double conj_BesK1 = gsl_sf_bessel_K1(conj_besk_inner_QX)/conj_X012;
    double upsilon_b_reg = Sq(z1)*(
        (2*z0*(z0 + z2) + Sq(z2))*(1 - 2*z1*(1-z1))*(conj_x_0p2sc1_dot_x_0p2sc1)*(conj_x02_dot_x02)/(conj_x02sq*x02sq)
        - z2 * (2*z0 + z2) * (2*z1 - 1) * (conj_x_0p2sc1_dot_conj_x20 * x_0p2sc1_dot_x20 - conj_x_0p2sc1_dot_x20 * x_0p2sc1_dot_conj_x20)/(conj_x02sq*x02sq)
    );
    double upsilon_c_reg = Sq(z0)*(
        (2*z1*(z1 + z2) + Sq(z2))*(1 - 2*z0*(1-z0))*(conj_x_0sc1p2_dot_x_0sc1p2)*(conj_x21_dot_x21)/(conj_x21sq*x21sq)
        - z2 * (2*z1 + z2) * (2*z0 - 1) * (conj_x_0sc1p2_dot_conj_x21 * x_0sc1p2_dot_x21 - conj_x_0sc1p2_dot_x21 * x_0sc1p2_dot_conj_x21)/(conj_x21sq*x21sq)
    );
    double upsilon_d = Sq(z0*z1*z2)/Sq(z0+z2) - Sq(z0*z1)*z1*z2/(z0+z2) * (x_0p2sc1_dot_x20/x02sq + conj_x_0p2sc1_dot_conj_x20/conj_x02sq)
                       + Sq(z0)*z1*Sq(z1+z2)*z2/(z0+z2) * (x_0sc1p2_dot_x21/x21sq + conj_x_0sc1p2_dot_conj_x21/conj_x21sq);
    double upsilon_e = Sq(z0*z1*z2)/Sq(z1+z2) - z0*Sq(z1*(z0+z2))*z2/(z1+z2) * (x_0p2sc1_dot_x20/x02sq + conj_x_0p2sc1_dot_conj_x20/conj_x02sq)
                       + z0*Sq(z0*z1)*z2/(z1+z2) * (x_0sc1p2_dot_x21/x21sq + conj_x_0sc1p2_dot_conj_x21/conj_x21sq);
    double upsilon_bc_interf = -z0*z1*( z1*(z0 + z2) + z0*(z1 + z2) ) *	( z0*(z0 + z2) + z1*(z1 + z2) )
                                * ( conj_x_0p2sc1_dot_x_0sc1p2 * (-1)*conj_x02_dot_x21/(conj_x02sq * x21sq) 
                                   + conj_x_0sc1p2_dot_x_0p2sc1 * (-1)*conj_x21_dot_x02 /(conj_x21sq * x02sq) )
                                + z0 * z1 * z2 * Sq(z0-z1) * (
                                    (conj_x_0p2sc1_dot_conj_x20 * x_0sc1p2_dot_x21 - conj_x_0p2sc1_dot_x21 * x_0sc1p2_dot_conj_x20)/(conj_x02sq * x21sq)
                                    +
                                    (conj_x_0sc1p2_dot_conj_x21 * x_0p2sc1_dot_x20 - conj_x_0sc1p2_dot_x20 * x_0p2sc1_dot_conj_x21)/(conj_x21sq * x02sq)
                                );
    double upsilon_terms = upsilon_b_reg + upsilon_c_reg + upsilon_d + upsilon_e + upsilon_bc_interf;
    double impact_fac_integrand = z0 * z1 * Sq(Q) * BesK1 * conj_BesK1 * upsilon_terms;
    
    res = I_Mx * impact_fac_integrand / (z0 * z1 * z2); // Remember to divide by the z_i under the integration measure
    return res;
}

