#include <gsl/gsl_sys.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include "nlodis_config.hpp"
#include "nlodissigmar_massiveq.hpp"


double ILdip_massive_LiLogConst(double Q, double z1, double x01sq, double mf) {
    // JP TODO
    double dip_res;
    return dip_res;
}

double ILdip_massive_Iab(double Q, double z1, double x01sq, double mf) {
    // JP TODO
    double dip_res;
    return dip_res;
}

double ILdip_massive_Icd(double Q, double z1, double x01sq, double mf) {
    // JP TODO
    double dip_res;
    return dip_res;
}

double ILNLOqg_massive(double Q, double x, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    // JP TODO
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double res = 0;
    if(gsl_finite(res)==1){
        return res;
    }else{
        res=0;
    }
    return res;
}