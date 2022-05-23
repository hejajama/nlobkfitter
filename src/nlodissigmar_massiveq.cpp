#include <gsl/gsl_sys.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_dilog.h>

#include "nlodis_config.hpp"
#include "nlodissigmar_massiveq.hpp"
#include "nlodissigmar.hpp"

// ----------------------- HELPER FUNCTIONS ------------------------

double L_dip( double Q, double z, double mf ) {
    // The L(gamma; z) function that appears in the longitudinal NLOdip part. 
    double gamma = sqrt( 1.0 + 4.0 * Sq(mf/Q) );

    double res = gsl_sf_dilog( 1.0 / ( 1.0 - 1.0 / (2.0 * z) * ( 1.0 - gamma) ) )
               + gsl_sf_dilog( 1.0 / ( 1.0 - 1.0 / (2.0 * z) * ( 1.0 + gamma) ) )
               + gsl_sf_dilog( 1.0 / ( 1.0 - 1.0 / (2.0 * (1.0-z)) * ( 1.0 - gamma) ) )
               + gsl_sf_dilog( 1.0 / ( 1.0 - 1.0 / (2.0 * (1.0-z)) * ( 1.0 + gamma) ) );

    return res;
}

double OmegaL_V( double Q, double z, double mf ) {
    // The Omega^L_V(gamma; z) function that appears in the longitudinal NLOdip part.
    double gamma = sqrt( 1.0 + 4.0 * Sq(mf/Q) );
    double res = 1.0/(2.0*z) * ( log( 1.0-z ) + gamma * log( (1.0+gamma)/(1.0+gamma-2.0*z) ) ) 
                +1.0/(2.0*(1.0-z)) * ( log( z ) + gamma * log( (1.0+gamma)/(1.0+gamma-2.0*(1.0-z)) ))
                +1.0/( 4.0*z*(1.0-z) ) * ( gamma-1.0 + 2.0 * Sq(mf/Q) ) * log( (z*(1.0-z)*Sq(Q) + Sq(mf) ) / Sq(mf) )  ;

    return res;
}


double OmegaT_V_unsymmetric( double Q, double z, double mf ) {
    // The first part of Omega^T_V(gamma; z) function that appears in the transverse NLOdip part.
    // Omega^T_V(gamma; z) = Omega^T_V_unsymmetric(gamma; z) + Omega^T_V_unsymmetric(gamma; 1-z)
    double gamma = sqrt( 1.0 + 4.0 * Sq(mf/Q) );
    double res = (1.0+1.0/(2.0*z)) * ( log(1.0-z) + gamma * log( (1.0+gamma)/(1.0+gamma-2.0*z) ) )
    - 1.0/(2.0*z) * ( (z + 0.5 ) * (1.0-gamma) + Sq(mf)/Sq(Q) ) * log( ( z*(1.0-z)*Sq(Q)+Sq(mf) )/Sq(mf) ) ;

    return res;
}


double OmegaT_N_unsymmetric( double Q, double z, double mf ) {
    // The first part of Omega^T_V(gamma; z) function that appears in the transverse NLOdip part.
    // Omega^T_N(gamma; z) = Omega^T_N_unsymmetric(gamma; z) - Omega^T_N_unsymmetric(gamma; 1-z)
    double gamma = sqrt( 1.0 + 4.0 * Sq(mf/Q) );
    double res = (1.0+z-2.0*Sq(z))/z * ( log(1.0-z) + gamma * log( (1.0+gamma)/(1.0+gamma-2.0*z) ) ) 
    + (1.0-z)/z * ( (0.5+z)*(gamma-1.0) - Sq(mf)/Sq(Q) ) * log(  (z*(1.0-z)*Sq(Q)+Sq(mf))/Sq(mf)) ;

    return res;
}


double IT_V1_unsymmetric( double Q, double z, double mf, double r, double xi ) {
    // The first part of the unintegrated I^T_V1 function that appears in the transverse NLOdip part.
    // I^T_V1 = I^T_V1_unsymmetric(z) + I^T_V1_unsymmetric(1-z)
    // Note that this has to be integrated over xi from 0 to 1.

    double kappa_z = sqrt( z*(1.0-z)*Sq(Q) + Sq(mf) );

    double term1 = 1.0/xi * ( 2.0*log(xi)/(1.0-xi) - (1.0+xi)/2.0 ) * ( sqrt( Sq(kappa_z) + xi/(1.0-xi) * (1.0-z) * Sq(mf)) * gsl_sf_bessel_K1( r*sqrt(Sq(kappa_z) + xi/(1.0-xi) *(1.0-z)* Sq(mf)) ) - kappa_z * gsl_sf_bessel_K1( r*kappa_z ) );
    double term2 = -( log(xi)/Sq(1.0-xi) + z/(1.0-xi) + z/2.0 ) * (1.0-z)*Sq(mf)/sqrt( Sq(kappa_z) + xi/(1.0-xi) * (1.0-z) *Sq(mf) ) * gsl_sf_bessel_K1( r*sqrt( Sq(kappa_z) + xi/(1.0-xi) * (1.0-z) *Sq(mf) ) );

    double res = term1 + term2;

    return res;

}


double IT_VMS1_unsymmetric( double Q, double z, double mf, double r, double xi ) {
    // The first part of the unintegrated I^T_VMS1 function that appears in the transverse NLOdip part.
    // I^T_VMS1 = I^T_VMS1_unsymmetric(z) + I^T_VMS1_unsymmetric(1-z)
    // Note that this has to be integrated over xi from 0 to 1.

    double kappa_z = sqrt( z*(1.0-z)*Sq(Q) + Sq(mf) );

    double term1 = 1.0/xi * (2.0 * log(xi)/(1.0-xi) - (1.0+xi)/(2.0)) * ( gsl_sf_bessel_K0( r*sqrt( Sq(kappa_z) + xi/(1.0-xi) * (1.0-z) *Sq(mf) ) ) - gsl_sf_bessel_K0(r*kappa_z) );
    double term2 = ( -3.0/2.0 * (1.0-z)/(1.0-xi) + (1.0-z)/2.0 ) * gsl_sf_bessel_K0( r* sqrt(Sq(kappa_z) + xi/(1.0-xi) * (1-z) *Sq(mf)) );

    double res = term1 + term2;

    return res;

}


double IT_V2_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u ) {
    // The first part of the unintegrated I^T_V2 function that appears in the transverse NLOdip part.
    // I^T_V2 = I^T_V2_unsymmetric(z) + I^T_V2_unsymmetric(1-z)
    // Note that this has to be integrated over y_chi and y_u, both from 0 to 1.

    double chi = z * y_chi;
    double u = (1.0-y_u)/y_u;

    double kappa_z = sqrt( z*(1.0-z)*Sq(Q) + Sq(mf) );
    double kappa_chi = sqrt( chi*(1.0-chi)*Sq(Q) + Sq(mf) );

    double term1 = - 1.0/(1.0-chi) * 1.0/(u*(u+1.0)) * Sq(mf)/Sq(kappa_chi) * ( 2.0*chi + Sq(  u/(u+1.0)) * 1.0/z * (z-chi) * (1.0-2.0*chi) ) * ( sqrt(Sq(kappa_z) + u*(1.0-z)/(1.0-chi)*Sq(kappa_chi)) * gsl_sf_bessel_K1(r*sqrt( Sq(kappa_z) + u*(1.0-z)/(1.0-chi)*Sq(kappa_chi) )) - kappa_z *gsl_sf_bessel_K1(r* kappa_z) );
    double term2 = -1.0/Sq(1.0-chi) * 1.0/(u+1.0) * (z-chi) * ( 1.0 - 2.0*u/(1.0+u)*(z-chi) + Sq(u/(u+1.0)) *1.0/z * Sq(z-chi) ) * Sq(mf)/sqrt( Sq(kappa_z) + u *(1.0-z)/(1.0-chi) * Sq(kappa_chi)) * gsl_sf_bessel_K1( r* sqrt( Sq(kappa_z) + u *(1.0-z)/(1.0-chi) * Sq(kappa_chi) ) );

    double jacobian = z / Sq(y_u);

    double res = jacobian * (term1 + term2);

    return res;


}
double IT_VMS2_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u ) {
    // The first part of the unintegrated I^T_VMS2 function that appears in the transverse NLOdip part.
    // I^T_VMS2 = I^T_VMS2_unsymmetric(z) + I^T_VMS2_unsymmetric(1-z)
    // Note that this has to be integrated over y_chi and y_u, both from 0 to 1.

    double chi = z * y_chi;
    double u = (1.0-y_u)/y_u;

    double kappa_z = sqrt( z*(1.0-z)*Sq(Q) + Sq(mf) );
    double kappa_chi = sqrt( chi*(1.0-chi)*Sq(Q) + Sq(mf) );

    double term1 = 1.0/(1.0-chi) * 1.0/Sq(u+1.0) * (-z - u/(1.0+u) * (z + u*chi)/z * (chi-(1.0-z))) * gsl_sf_bessel_K0(r * sqrt( Sq(kappa_z) + u*(1.0-z)/(1.0-chi) *Sq(kappa_chi) ));
    double term2 = 1.0/(u+1.0)/Sq(u+1.0) * ( Sq(kappa_z)/Sq(kappa_chi) * (1.0 + u * chi*(1.0-chi) / ( z*(1.0-z) )) - Sq(mf)/Sq(kappa_chi) * chi/(1.0-chi) * (2.0 *Sq(1.0+u)/u + u/(z*(1.0-z)) *Sq(z-chi) ) ) * ( gsl_sf_bessel_K0( r *sqrt( Sq(kappa_z) + u* (1.0-z)/(1.0-chi) *Sq(kappa_chi)) ) - gsl_sf_bessel_K0(r* kappa_z) );

    double jacobian = z / Sq(y_u);

    double res = jacobian * (term1 + term2);

    return res;

}


double IT_N_unsymmetric( double Q, double z, double mf, double r, double y_chi, double y_u ) {
    // The first part of the unintegrated I^T_N function that appears in the transverse NLOdip part.
    // I^T_N = I^T_N_unsymmetric(z) - I^T_N_unsymmetric(1-z)
    // Note that this has to be integrated over y_chi and y_u, both from 0 to 1.

    double chi = z * y_chi;
    double u = (1.0-y_u)/y_u;

    double kappa_z = sqrt( z*(1.0-z)*Sq(Q) + Sq(mf) );
    double kappa_chi = sqrt( chi*(1.0-chi)*Sq(Q) + Sq(mf) );

    double term1 = 2.0*(1.0-z)/z * 1.0/Sq(u+1.0)/(u+1.0) * ( (2.0+u)*u*z + Sq(u)*chi ) * sqrt( Sq(kappa_z) + u*(1.0-z)/(1.0-chi) * Sq(kappa_chi) ) * gsl_sf_bessel_K1( r* sqrt( Sq(kappa_z) + u*(1.0-z)/(1.0-chi) *Sq(kappa_chi) ) );
    double term2 = 2.0*(1.0-z)/z * 1.0/Sq(u+1.0)/(u+1.0) * Sq(mf)/Sq(kappa_chi) * ( z/(1.0-z) + chi/(1.0-chi) * (u-2.0*z-2.0*u*chi) ) * ( sqrt( Sq(kappa_z) + u*(1.0-z)/(1.0-chi)*Sq(kappa_chi)) * gsl_sf_bessel_K1( r * sqrt( Sq(kappa_z) + u*(1.0-z)/(1.0-chi) *Sq(kappa_chi) ) ) - kappa_z * gsl_sf_bessel_K1( r* kappa_z ) );

    double jacobian = z / Sq(y_u);

    double res = jacobian * (term1 + term2);

    return res;

}





struct G_qg_params{
    // Parameters for G_qg integration.
    int a, b;
    double Qbar, mf, x2, x3, omega, lambda;
};



int G_qg_integrand( const int *ndim, const double x[], const int *ncomp, double *f, void *params ) {
    // Integrand of G_qg.

    // Integration limits?????

    G_qg_params *p = (G_qg_params*) params;

    // Actual integration variables. 0 < y_i < 1
    double y_u = x[0];
    double y_t = x[1];


    int a = p->a;
    int b = p->b;
    double Qbar = p->Qbar;
    double mf = p->mf;
    double x2 = p->x2;
    double x3 = p->x3;
    double omega = p->omega;
    double lambda = p->lambda;

    double u = y_u / (1.0-y_u);
    double t, jacobian;

    if( omega > 1e-10 ) {
        double t_max = u / omega;
        t = y_t * t_max;
        jacobian = t_max / Sq(1.0-y_u) ; // Integration is performed over [0,1]^2 instead of [0,Inf]x[0,tmax].
    } else {
        t= y_t / (1.0-y_t);
        jacobian = 1.0 / Sq(1.0-y_u) / Sq(1.0-y_t);
    }

    double u_exp_inner_fun = u*(Sq(Qbar)+Sq(mf)) + Sq(x3)/(4.0*u);
    double t_exp_inner_fun = t*omega*lambda*Sq(mf) + Sq(x2)/(4.0*t);

    double u_part;
    double t_part;

    if( u_exp_inner_fun > 5e2 ) {
        u_part = 0.0;
    } else {
        u_part = 1.0/gsl_pow_int( u, a ) * gsl_sf_exp( -u_exp_inner_fun );
    }
    if( t_exp_inner_fun > 5e2 ) {
        t_part = 0.0;
    } else {
        t_part = 1.0/gsl_pow_int( t, b ) * gsl_sf_exp( -t_exp_inner_fun );
    }

    
    
    double res = u_part * t_part * jacobian;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}


double G_qg_int(int a, int b, double Qbar, double mf, double x2, double x3, double omega, double lambda) {
    // Functions G^(a;b)_(x) that appear in the qqg-part.

    string method = "cuhre";
    const int ndim = 2;


    double integral, error, prob;
    G_qg_params params;
    params.a = a;
    params.b = b;
    params.Qbar = Qbar;
    params.lambda = lambda;
    params.mf = mf;
    params.x2 = x2;
    params.x3 = x3;
    params.omega = omega;
    params.lambda = lambda;

    Cuba(method,ndim,G_qg_integrand,&params,&integral,&error,&prob);
    return integral;
}

double G_qg(int a, int b, double y_u, double y_t, double Qbar, double mf, double x2, double x3, double omega, double lambda) {
    // Functions G^(a;b)_(x) that appear in the qqg-part. Un-integrated form.

    double u = y_u / (1.0-y_u);
    double t, jacobian;

    if( omega > 1e-10 ) {
        double t_max = u / omega;
        t = y_t * t_max;
        jacobian = t_max / Sq(1.0-y_u) ; // Integration is performed over [0,1]^2 instead of [0,Inf]x[0,tmax].
    } else {
        t= y_t / (1.0-y_t);
        jacobian = 1.0 / Sq(1.0-y_u) / Sq(1.0-y_t);
    }

    double u_exp_inner_fun = u*(Sq(Qbar)+Sq(mf)) + Sq(x3)/(4.0*u);
    double t_exp_inner_fun = t*omega*lambda*Sq(mf) + Sq(x2)/(4.0*t);

    double u_part;
    double t_part;

    if( u_exp_inner_fun > 5e2 ) {
        u_part = 0.0;
    } else {
        u_part = 1.0/gsl_pow_int( u, a ) * gsl_sf_exp( -u_exp_inner_fun );
    }
    if( t_exp_inner_fun > 5e2 ) {
        t_part = 0.0;
    } else {
        t_part = 1.0/gsl_pow_int( t, b ) * gsl_sf_exp( -t_exp_inner_fun );
    }
    
    double res = u_part * t_part * jacobian;

    if(gsl_finite(res)==1){
        return res;
    }else{
        return 0;
    }
}


double IT_dipole_jk_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    // double omega_j = z0*z2/(z1*Sq(z0+z2));
    // double omega_k = z1*z2/(z0*Sq(z1+z2));
    // double lambda_j = z1*z2/z0;
    // double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double term_j = 1.0/Sq(z0+z2) * (2.0*z0*(z0+z2)+Sq(z2)) * ( 1.0-2.0*z1*(1.0-z1) ) * ( Sq(Qbar_j) + Sq(mf) ) / Sq(x2_j) 
                    * ( - exp( -Sq(x2_j) / ( x01sq *exp(M_EULER) ) ) * Sq(gsl_sf_bessel_K1( sqrt( x01sq * ( Sq(Qbar_j) + Sq(mf) ) ) ) ));
    double term_k = 1.0/Sq(z1+z2) * (2.0*z1*(z1+z2)+Sq(z2)) * ( 1.0-2.0*z0*(1.0-z0) ) * ( Sq(Qbar_k) + Sq(mf) ) / Sq(x2_k) 
                    * ( - exp( -Sq(x2_k) / ( x01sq *exp(M_EULER) ) ) * Sq(gsl_sf_bessel_K1( sqrt( x01sq * ( Sq(Qbar_k) + Sq(mf) ) ) ) ));
   
    double res = term_j + term_k;

    return res;
}



double IT_tripole_jk_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    // double lambda_j = z1*z2/z0;
    // double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double term_j = 1.0/Sq(z0+z2) * (2.0*z0*(z0+z2)+Sq(z2)) * ( 1.0-2.0*z1*(1.0-z1) ) * ( Sq(Qbar_j) + Sq(mf) ) / Sq(x2_j) 
                    * Sq(x3_j) / ( Sq(x3_j) + omega_j * Sq(x2_j) ) * Sq( gsl_sf_bessel_K1( sqrt( Sq(x3_j) + omega_j * Sq(x2_j) ) * sqrt( Sq(Qbar_j) + Sq(mf) ) ) ) ;
    double term_k = 1.0/Sq(z1+z2) * (2.0*z1*(z1+z2)+Sq(z2)) * ( 1.0-2.0*z0*(1.0-z0) ) * ( Sq(Qbar_k) + Sq(mf) ) / Sq(x2_k) 
                    * Sq(x3_k) / ( Sq(x3_k) + omega_k * Sq(x2_k) ) * Sq( gsl_sf_bessel_K1( sqrt( Sq(x3_k) + omega_k * Sq(x2_k) ) * sqrt( Sq(Qbar_k) + Sq(mf) ) ) ) ;
   
    double res = term_j + term_k;

    return res;
}


double IT_tripole_jk_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;

    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double u = (1-y_u)/y_u;
    double t_j = y_t * u / omega_j;
    double t_k = y_t * u / omega_k;

    double g_bar_j = exp( -u * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u) ) * exp( -x02sq/(4.0*t_j) ) * ( exp(-t_j*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k = exp( -u * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u) ) * exp( -x21sq/(4.0*t_k) ) * ( exp(-t_k*omega_k*lambda_k*Sq(mf)) - 1.0 );


    double term_j = 1.0/(Sq(y_u) * u * omega_j * Sq(t_j)  )
                    * 1.0/Sq(z0+z2) * (2.0*z0*(z0+z2)+Sq(z2)) * ( 1.0-2.0*z1*(1.0-z1) ) 
                    * g_bar_j * Sq(x3_j)/8.0 * sqrt( ( Sq(Qbar_j) + Sq(mf) ) / (Sq(x3_j) + omega_j * Sq(x2_j) ) )
                    * gsl_sf_bessel_K1( sqrt( Sq(x3_j) + omega_j * Sq(x2_j) ) * sqrt( Sq(Qbar_j) + Sq(mf) ) )  ;
    double term_k = 1.0/(Sq(y_u) * u * omega_k * Sq(t_k)  )
                    * 1.0/Sq(z1+z2) * (2.0*z1*(z1+z2)+Sq(z2)) * ( 1.0-2.0*z0*(1.0-z0) ) 
                    * g_bar_k * Sq(x3_k)/8.0 * sqrt( ( Sq(Qbar_k) + Sq(mf) ) / (Sq(x3_k) + omega_k * Sq(x2_k) ) )
                    * gsl_sf_bessel_K1( sqrt( Sq(x3_k) + omega_k * Sq(x2_k) ) * sqrt( Sq(Qbar_k) + Sq(mf) ) )  ;

    double res = term_j + term_k;

    return res;
}




double IT_tripole_jk_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;

    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double u1 = (1-y_u1)/y_u1;
    double u2 = (1-y_u2)/y_u2;
    double t_j1 = y_t1 * u1 / omega_j;
    double t_k1 = y_t1 * u1 / omega_k;
    double t_j2 = y_t2 * u2 / omega_j;
    double t_k2 = y_t2 * u2 / omega_k;

    double g_bar_j1 = exp( -u1 * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u1) ) * exp( -x02sq/(4.0*t_j1) ) * ( exp(-t_j1*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k1 = exp( -u1 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u1) ) * exp( -x21sq/(4.0*t_k1) ) * ( exp(-t_k1*omega_k*lambda_k*Sq(mf)) - 1.0 );
    double g_bar_j2 = exp( -u2 * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u2) ) * exp( -x02sq/(4.0*t_j2) ) * ( exp(-t_j2*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k2 = exp( -u2 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u2) ) * exp( -x21sq/(4.0*t_k2) ) * ( exp(-t_k2*omega_k*lambda_k*Sq(mf)) - 1.0 );


    double term_j = 1.0/(Sq(y_u1) * u1 * omega_j * Sq(t_j1)  ) * 1.0/(Sq(y_u2) * u2 * omega_j * Sq(t_j2)  )
                    * 1.0/Sq(z0+z2) * (2.0*z0*(z0+z2)+Sq(z2)) * ( 1.0-2.0*z1*(1.0-z1) ) 
                    * Sq(x3_j) * Sq(x2_j) / 256.0 * g_bar_j1 * g_bar_j2 ;
    double term_k = 1.0/(Sq(y_u1) * u1 * omega_k * Sq(t_k1)  ) * 1.0/(Sq(y_u2) * u2 * omega_k * Sq(t_k2)  )
                    * 1.0/Sq(z1+z2) * (2.0*z1*(z1+z2)+Sq(z2)) * ( 1.0-2.0*z0*(1.0-z0) ) 
                    * Sq(x3_k) * Sq(x2_k) / 256.0 * g_bar_k1 * g_bar_k2 ;

    double res = term_j + term_k;

    return res;
}


double IT_dipole_jkm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    // double omega_j = z0*z2/(z1*Sq(z0+z2));
    // double omega_k = z1*z2/(z0*Sq(z1+z2));
    // double lambda_j = z1*z2/z0;
    // double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double term_j = 1.0/Sq(z0+z2) * (2.0*z0*(z0+z2)+Sq(z2)) / Sq(x2_j) 
                    * ( - exp( -Sq(x2_j) / ( x01sq *exp(M_EULER) ) ) * Sq(gsl_sf_bessel_K0( sqrt( x01sq * ( Sq(Qbar_j) + Sq(mf) ) ) ) ));
    double term_k = 1.0/Sq(z1+z2) * (2.0*z1*(z1+z2)+Sq(z2)) / Sq(x2_k) 
                    * ( - exp( -Sq(x2_k) / ( x01sq *exp(M_EULER) ) ) * Sq(gsl_sf_bessel_K0( sqrt( x01sq * ( Sq(Qbar_k) + Sq(mf) ) ) ) ));
   
    double res = Sq(mf) * (term_j + term_k);

    return res;
}

double IT_tripole_jkm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    // double lambda_j = z1*z2/z0;
    // double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double term_j = 1.0/Sq(z0+z2) * (2.0*z0*(z0+z2)+Sq(z2)) / Sq(x2_j) 
                    * Sq( gsl_sf_bessel_K0( sqrt( Sq(x3_j) + omega_j * Sq(x2_j) ) * sqrt( Sq(Qbar_j) + Sq(mf) ) ) );
    double term_k = 1.0/Sq(z1+z2) * (2.0*z1*(z1+z2)+Sq(z2)) / Sq(x2_k) 
                    * Sq( gsl_sf_bessel_K0( sqrt( Sq(x3_k) + omega_k * Sq(x2_k) ) * sqrt( Sq(Qbar_k) + Sq(mf) ) ) );
   
    double res = Sq(mf) * (term_j + term_k);

    return res;
}



double IT_tripole_jkm_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;

    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double u = (1-y_u)/y_u;
    double t_j = y_t * u / omega_j;
    double t_k = y_t * u / omega_k;

    double g_bar_j = exp( -u * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u) ) * exp( -x02sq/(4.0*t_j) ) * ( exp(-t_j*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k = exp( -u * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u) ) * exp( -x21sq/(4.0*t_k) ) * ( exp(-t_k*omega_k*lambda_k*Sq(mf)) - 1.0 );


    double term_j = 1.0/(Sq(y_u) * omega_j * Sq(t_j)  )
                    * 1.0/Sq(z0+z2) * (2.0*z0*(z0+z2)+Sq(z2))
                    * g_bar_j/4.0 * gsl_sf_bessel_K0( sqrt( Sq(x3_j) + omega_j * Sq(x2_j) ) * sqrt( Sq(Qbar_j) + Sq(mf) ) )  ;
    double term_k = 1.0/(Sq(y_u) * omega_k * Sq(t_k)  )
                    * 1.0/Sq(z1+z2) * (2.0*z1*(z1+z2)+Sq(z2))
                    * g_bar_k/4.0 * gsl_sf_bessel_K0( sqrt( Sq(x3_k) + omega_k * Sq(x2_k) ) * sqrt( Sq(Qbar_k) + Sq(mf) ) )  ;

    double res = Sq(mf) * (term_j + term_k);

    return res;
}


double IT_tripole_jkm_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;

    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double u1 = (1-y_u1)/y_u1;
    double u2 = (1-y_u2)/y_u2;
    double t_j1 = y_t1 * u1 / omega_j;
    double t_k1 = y_t1 * u1 / omega_k;
    double t_j2 = y_t2 * u2 / omega_j;
    double t_k2 = y_t2 * u2 / omega_k;

    double g_bar_j1 = exp( -u1 * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u1) ) * exp( -x02sq/(4.0*t_j1) ) * ( exp(-t_j1*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k1 = exp( -u1 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u1) ) * exp( -x21sq/(4.0*t_k1) ) * ( exp(-t_k1*omega_k*lambda_k*Sq(mf)) - 1.0 );
    double g_bar_j2 = exp( -u2 * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u2) ) * exp( -x02sq/(4.0*t_j2) ) * ( exp(-t_j2*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k2 = exp( -u2 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u2) ) * exp( -x21sq/(4.0*t_k2) ) * ( exp(-t_k2*omega_k*lambda_k*Sq(mf)) - 1.0 );


    double term_j = 1.0/(Sq(y_u1) *  omega_j * Sq(t_j1)  ) * 1.0/(Sq(y_u2) *  omega_j * Sq(t_j2)  )
                    * 1.0/Sq(z0+z2) * (2.0*z0*(z0+z2)+Sq(z2)) 
                    * Sq(x2_j) / 64.0 * g_bar_j1 * g_bar_j2 ;
    double term_k = 1.0/(Sq(y_u1) *  omega_k * Sq(t_k1)  ) * 1.0/(Sq(y_u2) *  omega_k * Sq(t_k2)  )
                    * 1.0/Sq(z1+z2) * (2.0*z1*(z1+z2)+Sq(z2)) 
                    * Sq(x2_k) / 64.0 * g_bar_k1 * g_bar_k2 ;

    double res = Sq(mf) * (term_j + term_k);

    return res;
}




double IT_tripole_F_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double x2j_x3j = x20x21 - z0/(z0+z2) * x02sq ;
    double x2k_x3k = -x20x21 + z1/(z1+z2) * x21sq ;
    double x2j_x3k = -x02sq + z1/(z1+z2) * x20x21 ;
    double x2k_x3j = x21sq - z0/(z0+z2) * x20x21 ;
    double x3j_x3k = z0/(z0+z2)*x02sq + z1/(z1+z2)*x21sq - ( 1.0 + z0*z1/(z0+z2)/(z1+z2) ) * x20x21;

    double G22_sing_j = 1.0/Sq(x2_j) * sqrt( ( Sq(Qbar_j) + Sq(mf) ) / ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( Sq(Qbar_j) + Sq(mf) ) * ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) );
    double G22_sing_k = 1.0/Sq(x2_k) * sqrt( ( Sq(Qbar_k) + Sq(mf) ) / ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( Sq(Qbar_k) + Sq(mf) ) * ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) );

    double H_j = 4.0 * sqrt( ( Sq(Qbar_j) + Sq(mf) * (1.0+lambda_j) ) / ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( Sq(Qbar_j) + Sq(mf) * (1.0+lambda_j) ) * ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) );
    double H_k = 4.0 * sqrt( ( Sq(Qbar_k) + Sq(mf) * (1.0+lambda_k) ) / ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( Sq(Qbar_k) + Sq(mf) * (1.0+lambda_k) ) * ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) );

    double term_1 = 4.0/(z0+z2)/(z1+z2) * ( z2*Sq(z0-z1) * ( x2j_x3j * x2k_x3k - x2k_x3j * x2j_x3k ) 
                    - ( z1*(z0+z2) + z0*(z1+z2) )*(z0*(z0+z2)+z1*(z1+z2)) * x20x21 * x3j_x3k ) *
                    G22_sing_j * G22_sing_k;
    double term_2j = -(z0+z2)*z1*z2/Sq(z1+z2) * x2j_x3j * H_k * G22_sing_j;
    double term_2k = (z1+z2)*z0*z2/Sq(z0+z2) * x2k_x3k * H_j * G22_sing_k;
    double term_3j = -Sq(z0)*z1*z2/(z0+z2)/Sq(z0+z2) * x2j_x3j * H_j * G22_sing_j ;
    double term_3k = Sq(z1)*z0*z2/(z1+z2)/Sq(z1+z2) * x2k_x3k * H_k * G22_sing_k ;
    double term_4j = Sq(z0*z2)/(8.0*Sq(z0+z2)*Sq(z0+z2))*Sq(H_j);
    double term_4k = Sq(z1*z2)/(8.0*Sq(z1+z2)*Sq(z1+z2))*Sq(H_k);

    double res = 0.5 * (term_1 + term_2j + term_2k + term_3j + term_3k + term_4j + term_4k);

    return res;
}


double IT_tripole_F_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double x2j_x3j = x20x21 - z0/(z0+z2) * x02sq ;
    double x2k_x3k = -x20x21 + z1/(z1+z2) * x21sq ;
    double x2j_x3k = -x02sq + z1/(z1+z2) * x20x21 ;
    double x2k_x3j = x21sq - z0/(z0+z2) * x20x21 ;
    double x3j_x3k = z0/(z0+z2)*x02sq + z1/(z1+z2)*x21sq - ( 1.0 + z0*z1/(z0+z2)/(z1+z2) ) * x20x21;


    double u = (1-y_u)/y_u;
    double t_j = y_t * u / omega_j;
    double t_k = y_t * u / omega_k;

    double g_bar_j = exp( -u * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u) ) * exp( -x02sq/(4.0*t_j) ) * ( exp(-t_j*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k = exp( -u * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u) ) * exp( -x21sq/(4.0*t_k) ) * ( exp(-t_k*omega_k*lambda_k*Sq(mf)) - 1.0 );


    double G22_sing_j = 1.0/Sq(x2_j) * sqrt( ( Sq(Qbar_j) + Sq(mf) ) / ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( Sq(Qbar_j) + Sq(mf) ) * ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) );
    double G22_sing_k = 1.0/Sq(x2_k) * sqrt( ( Sq(Qbar_k) + Sq(mf) ) / ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( Sq(Qbar_k) + Sq(mf) ) * ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) );

    double H_j = 4.0 * sqrt( ( Sq(Qbar_j) + Sq(mf) * (1.0+lambda_j) ) / ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( Sq(Qbar_j) + Sq(mf) * (1.0+lambda_j) ) * ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) );
    double H_k = 4.0 * sqrt( ( Sq(Qbar_k) + Sq(mf) * (1.0+lambda_k) ) / ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( Sq(Qbar_k) + Sq(mf) * (1.0+lambda_k) ) * ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) );

    double term_1 = 1.0/(4.0*(z0+z2)*(z1+z2)) * (z2*Sq(z0-z1)* ( x2j_x3j * x2k_x3k - x2k_x3j * x2j_x3k ) 
                    - (z1*(z0+z2)+z0*(z1+z2))* (z0*(z0+z2)+z1*(z1+z2))*x20x21 * x3j_x3k  )
                    * (
                        1.0/(Sq(y_u)*u*omega_k*Sq(t_k))*g_bar_k * G22_sing_j
                        +1.0/(Sq(y_u)*u*omega_j*Sq(t_j))*g_bar_j * G22_sing_k
                    );
    double term_2j = -(z0+z2)*z1*z2/(16.0*Sq(z1+z2)) * x2j_x3j * H_k *1.0/(Sq(y_u)*u*omega_j*Sq(t_j))*g_bar_j;
    double term_2k = (z1+z2)*z0*z2/(16.0*Sq(z0+z2)) * x2k_x3k * H_j *1.0/(Sq(y_u)*u*omega_k*Sq(t_k))*g_bar_k;
    double term_3j = -Sq(z0)*z1*z2/(16.0*(z0+z2)*Sq(z0+z2)) * x2j_x3j*H_j*1.0/(Sq(y_u)*u*omega_j*Sq(t_j))*g_bar_j;
    double term_3k = Sq(z1)*z0*z2/(16.0*(z1+z2)*Sq(z1+z2)) * x2k_x3k*H_k*1.0/(Sq(y_u)*u*omega_k*Sq(t_k))*g_bar_k;

    double res = 0.5 * (term_1 + term_2j + term_2k + term_3j + term_3k);

    return res;
}


double IT_tripole_F_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;

    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double x2j_x3j = x20x21 - z0/(z0+z2) * x02sq ;
    double x2k_x3k = -x20x21 + z1/(z1+z2) * x21sq ;
    double x2j_x3k = -x02sq + z1/(z1+z2) * x20x21 ;
    double x2k_x3j = x21sq - z0/(z0+z2) * x20x21 ;
    double x3j_x3k = z0/(z0+z2)*x02sq + z1/(z1+z2)*x21sq - ( 1.0 + z0*z1/(z0+z2)/(z1+z2) ) * x20x21;

    double u1 = (1-y_u1)/y_u1;
    double u2 = (1-y_u2)/y_u2;
    double t_j1 = y_t1 * u1 / omega_j;
    double t_k1 = y_t1 * u1 / omega_k;
    double t_j2 = y_t2 * u2 / omega_j;
    double t_k2 = y_t2 * u2 / omega_k;

    double g_bar_j1 = exp( -u1 * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u1) ) * exp( -x02sq/(4.0*t_j1) ) * ( exp(-t_j1*omega_j*lambda_j*Sq(mf)) - 1.0 );
    // double g_bar_k1 = exp( -u1 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u1) ) * exp( -x21sq/(4.0*t_k1) ) * ( exp(-t_k1*omega_k*lambda_k*Sq(mf)) - 1.0 );
    // double g_bar_j2 = exp( -u2 * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u2) ) * exp( -x02sq/(4.0*t_j2) ) * ( exp(-t_j2*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k2 = exp( -u2 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u2) ) * exp( -x21sq/(4.0*t_k2) ) * ( exp(-t_k2*omega_k*lambda_k*Sq(mf)) - 1.0 );

    double res = 0.5 * 1.0/( Sq(y_u1) * u1 * omega_j * Sq(t_j1) * Sq(y_u2) * u2 * omega_k * Sq(t_k2) ) * g_bar_j1 * g_bar_k2
                 *1.0/(64.0*(z0+z2)*(z1+z2)) * (z2*Sq(z0-z1)* ( x2j_x3j * x2k_x3k - x2k_x3j * x2j_x3k ) 
                    - (z1*(z0+z2)+z0*(z1+z2))* (z0*(z0+z2)+z1*(z1+z2))*x20x21 * x3j_x3k  );

    return res;
}



double IT_tripole_Fm_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    // double x2j_x3j = x20x21 - z0/(z0+z2) * x02sq ;
    // double x2k_x3k = -x20x21 + z1/(z1+z2) * x21sq ;
    // double x2j_x3k = -x02sq + z1/(z1+z2) * x20x21 ;
    // double x2k_x3j = x21sq - z0/(z0+z2) * x20x21 ;
    // double x3j_x3k = z0/(z0+z2)*x02sq + z1/(z1+z2)*x21sq - ( 1.0 + z0*z1/(z0+z2)/(z1+z2) ) * x20x21;

    double G12_sing_j = 1.0/Sq(x2_j) * gsl_sf_bessel_K0( sqrt(( Sq(Qbar_j) + Sq(mf) ) * ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) );
    double G12_sing_k = 1.0/Sq(x2_k) * gsl_sf_bessel_K0( sqrt(( Sq(Qbar_k) + Sq(mf) ) * ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) );

    // double G22_sing_j = 1.0/Sq(x2_j) * sqrt( ( Sq(Qbar_j) + Sq(mf) ) / ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) 
    //                     * gsl_sf_bessel_K1( sqrt(( Sq(Qbar_j) + Sq(mf) ) * ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) );
    // double G22_sing_k = 1.0/Sq(x2_k) * sqrt( ( Sq(Qbar_k) + Sq(mf) ) / ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) 
    //                     * gsl_sf_bessel_K1( sqrt(( Sq(Qbar_k) + Sq(mf) ) * ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) );

    // double H_j = 4.0 * sqrt( ( Sq(Qbar_j) + Sq(mf) * (1.0+lambda_j) ) / ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) 
    //             * gsl_sf_bessel_K1( sqrt( ( Sq(Qbar_j) + Sq(mf) * (1.0+lambda_j) ) * ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) );
    // double H_k = 4.0 * sqrt( ( Sq(Qbar_k) + Sq(mf) * (1.0+lambda_k) ) / ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) 
    //             * gsl_sf_bessel_K1( sqrt( ( Sq(Qbar_k) + Sq(mf) * (1.0+lambda_k) ) * ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) );



    double res = 0.5*Sq(mf) * (
        -1.0/(32.0*(z0+z2)*(z1+z2)) * ( (2.0*z0+z2)*(2.0*z1+z2)+Sq(z2) ) *x20x21 * 8.0 * G12_sing_j * 8.0 *G12_sing_k
    );

    return res;
}


double IT_tripole_Fm_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;


    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double x2j_x3j = x20x21 - z0/(z0+z2) * x02sq ;
    double x2k_x3k = -x20x21 + z1/(z1+z2) * x21sq ;
    double x2j_x3k = -x02sq + z1/(z1+z2) * x20x21 ;
    double x2k_x3j = x21sq - z0/(z0+z2) * x20x21 ;
    double x3j_x3k = z0/(z0+z2)*x02sq + z1/(z1+z2)*x21sq - ( 1.0 + z0*z1/(z0+z2)/(z1+z2) ) * x20x21;


    double u = (1-y_u)/y_u;
    double t_j = y_t * u / omega_j;
    double t_k = y_t * u / omega_k;

    double g_bar_j = exp( -u * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u) ) * exp( -x02sq/(4.0*t_j) ) * ( exp(-t_j*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k = exp( -u * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u) ) * exp( -x21sq/(4.0*t_k) ) * ( exp(-t_k*omega_k*lambda_k*Sq(mf)) - 1.0 );


    double g_j = exp( -u * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u) ) * exp( -x02sq/(4.0*t_j) ) * exp(-t_j*omega_j*lambda_j*Sq(mf));
    double g_k = exp( -u * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u) ) * exp( -x21sq/(4.0*t_k) ) * exp(-t_k*omega_k*lambda_k*Sq(mf));

    double G12_sing_j = 1.0/Sq(x2_j) * gsl_sf_bessel_K0( sqrt(( Sq(Qbar_j) + Sq(mf) ) * ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) );
    double G12_sing_k = 1.0/Sq(x2_k) * gsl_sf_bessel_K0( sqrt(( Sq(Qbar_k) + Sq(mf) ) * ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) );

    double G22_sing_j = 1.0/Sq(x2_j) * sqrt( ( Sq(Qbar_j) + Sq(mf) ) / ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( Sq(Qbar_j) + Sq(mf) ) * ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) );
    double G22_sing_k = 1.0/Sq(x2_k) * sqrt( ( Sq(Qbar_k) + Sq(mf) ) / ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) 
                        * gsl_sf_bessel_K1( sqrt(( Sq(Qbar_k) + Sq(mf) ) * ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) );

    double H_j = 4.0 * sqrt( ( Sq(Qbar_j) + Sq(mf) * (1.0+lambda_j) ) / ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( Sq(Qbar_j) + Sq(mf) * (1.0+lambda_j) ) * ( Sq(x3_j) + omega_j * Sq(x2_j) ) ) );
    double H_k = 4.0 * sqrt( ( Sq(Qbar_k) + Sq(mf) * (1.0+lambda_k) ) / ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) 
                * gsl_sf_bessel_K1( sqrt( ( Sq(Qbar_k) + Sq(mf) * (1.0+lambda_k) ) * ( Sq(x3_k) + omega_k * Sq(x2_k) ) ) );

    double term_1j = -z0*z1*Sq(z2)/(16.0*(z0+z2)*Sq(z0+z2)) * x2j_x3j * u/(Sq(y_u)*omega_j*Sq(u)*t_j) * g_j * 8.0 * G12_sing_j;
    double term_1k =  z0*z1*Sq(z2)/(16.0*(z1+z2)*Sq(z1+z2)) * x2k_x3k * u/(Sq(y_u)*omega_k*Sq(u)*t_k) * g_k * 8.0 * G12_sing_k;
    double term_2 = -1.0/(32.0*(z0+z2)*(z1+z2))* ( (2.0*z0+z2)*(2.0*z1+z2) + Sq(z2) ) * x20x21 * (
         u/(Sq(y_u)*omega_k*u*Sq(t_k)) * g_bar_k * 8.0 * G12_sing_j
        +u/(Sq(y_u)*omega_j*u*Sq(t_j)) * g_bar_j * 8.0 * G12_sing_k
    );
    double term_3j = -Sq(z0*z2)/(16.0*(z0+z2)*Sq(z1+z2)) * x2j_x3k * u/(Sq(y_u)*omega_k *Sq(u)*t_k ) * g_k * 8.0 * G12_sing_j;
    double term_3k =  Sq(z1*z2)/(16.0*Sq(z0+z2)*(z1+z2)) * x2k_x3j * u/(Sq(y_u)*omega_j *Sq(u)*t_j ) * g_j * 8.0 * G12_sing_k;
    double term_4j = -z0*z1*Sq(z2)/(16.0*(z0+z2)*Sq(z0+z2)) * x2j_x3j * u/(Sq(y_u)*omega_j *u*t_j ) * g_j * 16.0 * G22_sing_j ;
    double term_4k =  z0*z1*Sq(z2)/(16.0*(z1+z2)*Sq(z1+z2)) * x2k_x3k * u/(Sq(y_u)*omega_k *u*t_k ) * g_k * 16.0 * G22_sing_k ;
    double term_5j = -(z0+z2)*Sq(z2)/(16.0*Sq(z1+z2)) * x2j_x3j * u/(Sq(y_u)*omega_k *u*t_k ) * g_k * 16.0 * G22_sing_j ;
    double term_5k =  (z1+z2)*Sq(z2)/(16.0*Sq(z0+z2)) * x2k_x3k * u/(Sq(y_u)*omega_j *u*t_j ) * g_j * 16.0 * G22_sing_k ;
    double term_6j = z0*z2*Sq(z2)/(4.0*Sq(z0+z2)*Sq(z0+z2))* H_j * u/(Sq(y_u)*omega_j *u*t_j ) * g_j;
    double term_6k = z1*z2*Sq(z2)/(4.0*Sq(z1+z2)*Sq(z1+z2))* H_k * u/(Sq(y_u)*omega_k *u*t_k ) * g_k;

    double res = 0.5*Sq(mf) * ( term_1j + term_1k + term_2 + term_3j + term_3k + term_4j + term_4k + term_5j + term_5k + term_6j + term_6k );

    return res;
}

double IT_tripole_Fm_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2){

    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_j = Q*sqrt(z1*(1.0-z1));
    double Qbar_k = Q*sqrt(z0*(1.0-z0));
    double omega_j = z0*z2/(z1*Sq(z0+z2));
    double omega_k = z1*z2/(z0*Sq(z1+z2));
    double lambda_j = z1*z2/z0;
    double lambda_k = z0*z2/z1;

    double x2_j = sqrt(x02sq);
    double x2_k = sqrt(x21sq);
    double x3_j = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_k = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double x2j_x3j = x20x21 - z0/(z0+z2) * x02sq ;
    double x2k_x3k = -x20x21 + z1/(z1+z2) * x21sq ;
    double x2j_x3k = -x02sq + z1/(z1+z2) * x20x21 ;
    double x2k_x3j = x21sq - z0/(z0+z2) * x20x21 ;
    double x3j_x3k = z0/(z0+z2)*x02sq + z1/(z1+z2)*x21sq - ( 1.0 + z0*z1/(z0+z2)/(z1+z2) ) * x20x21;

    double u1 = (1-y_u1)/y_u1;
    double u2 = (1-y_u2)/y_u2;
    double t_j1 = y_t1 * u1 / omega_j;
    double t_k1 = y_t1 * u1 / omega_k;
    double t_j2 = y_t2 * u2 / omega_j;
    double t_k2 = y_t2 * u2 / omega_k;

    double g_bar_j1 = exp( -u1 * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u1) ) * exp( -x02sq/(4.0*t_j1) ) * ( exp(-t_j1*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k1 = exp( -u1 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u1) ) * exp( -x21sq/(4.0*t_k1) ) * ( exp(-t_k1*omega_k*lambda_k*Sq(mf)) - 1.0 );
    double g_bar_j2 = exp( -u2 * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u2) ) * exp( -x02sq/(4.0*t_j2) ) * ( exp(-t_j2*omega_j*lambda_j*Sq(mf)) - 1.0 );
    double g_bar_k2 = exp( -u2 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u2) ) * exp( -x21sq/(4.0*t_k2) ) * ( exp(-t_k2*omega_k*lambda_k*Sq(mf)) - 1.0 );

    double g_j1 = exp( -u1 * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u1) ) * exp( -x02sq/(4.0*t_j1) ) * exp(-t_j1*omega_j*lambda_j*Sq(mf));
    double g_k1 = exp( -u1 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u1) ) * exp( -x21sq/(4.0*t_k1) ) * exp(-t_k1*omega_k*lambda_k*Sq(mf));
    double g_j2 = exp( -u2 * (Sq(Qbar_j)+Sq(mf)) - Sq(x3_j)/(4.0*u2) ) * exp( -x02sq/(4.0*t_j2) ) * exp(-t_j2*omega_j*lambda_j*Sq(mf));
    double g_k2 = exp( -u2 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u2) ) * exp( -x21sq/(4.0*t_k2) ) * exp(-t_k2*omega_k*lambda_k*Sq(mf));

    double term_1j = Sq(z2)*Sq(z2)/(64.0*Sq(z0+z2)*Sq(z0+z2)) * ( 4.0*z1*(z1-1.0)+2.0 )*Sq(x3_j) 
                    * u1/( Sq(y_u1)*omega_j *Sq(u1)*t_j1) *u2/( Sq(y_u2)*omega_j *Sq(u2)*t_j2) * g_j1 * g_j2 ;
    double term_1k = Sq(z2)*Sq(z2)/(64.0*Sq(z1+z2)*Sq(z1+z2)) * ( 4.0*z0*(z0-1.0)+2.0 )*Sq(x3_k) 
                    * u1/( Sq(y_u1)*omega_k *Sq(u1)*t_k1) *u2/( Sq(y_u2)*omega_k *Sq(u2)*t_k2) * g_k1 * g_k2 ;
    double term_2j = -z0*z1*Sq(z2)/(16.0*(z0+z2)*Sq(z0+z2))* x2j_x3j
                    * u1/( Sq(y_u1)*omega_j *u1*Sq(t_j1)) *u2/( Sq(y_u2)*omega_j *Sq(u2)*t_j2) * g_bar_j1 * g_j2 ;
    double term_2k =  z0*z1*Sq(z2)/(16.0*(z1+z2)*Sq(z1+z2))* x2k_x3k
                    * u1/( Sq(y_u1)*omega_k *u1*Sq(t_k1)) *u2/( Sq(y_u2)*omega_k *Sq(u2)*t_k2) * g_bar_k1 * g_k2 ;
    double term_3a = -1.0/(32.0*(z0+z2)*(z1+z2)) * ( (2.0*z0+z2)*(2.0*z1+z2) +Sq(z2) ) * x20x21
                    * u1/( Sq(y_u1)*omega_j *u1*Sq(t_j1)) *u2/( Sq(y_u2)*omega_k *u2*Sq(t_k2)) * g_bar_j1 * g_bar_k2 ;
    double term_3b = Sq(z2)*Sq(z2)/(32.0*Sq(z0+z2)*Sq(z1+z2)) * ( (2.0*z0+z2)*(2.0*z1+z2) + Sq(z2) ) * x3j_x3k
                    * u1/( Sq(y_u1)*omega_j *Sq(u1)*t_j1) *u2/( Sq(y_u2)*omega_k *Sq(u2)*t_k2) * g_j1 * g_k2 ;
    double term_4j = Sq(mf)/16.0 * 2.0*Sq(z2/(z0+z2))*Sq(z2/(z0+z2)) * u1/( Sq(y_u1)*omega_j *u1*t_j1) *u2/( Sq(y_u2)*omega_j *u2*t_j2) * g_j1 * g_j2 ;
    double term_4k = Sq(mf)/16.0 * 2.0*Sq(z2/(z1+z2))*Sq(z2/(z1+z2)) * u1/( Sq(y_u1)*omega_k *u1*t_k1) *u2/( Sq(y_u2)*omega_k *u2*t_k2) * g_k1 * g_k2 ;
    double term_5j = -Sq(z0*z2)/(16.0*(z0+z2)*Sq(z1+z2))*x2j_x3k
                    * u1/( Sq(y_u1)*omega_j *u1*Sq(t_j1)) *u2/( Sq(y_u2)*omega_k *Sq(u2)*t_k2) * g_bar_j1 * g_k2 ;
    double term_5k =  Sq(z1*z2)/(16.0*(z1+z2)*Sq(z0+z2))*x2k_x3j
                    * u1/( Sq(y_u1)*omega_k *u1*Sq(t_k1)) *u2/( Sq(y_u2)*omega_j *Sq(u2)*t_j2) * g_bar_k1 * g_j2 ;
    double term_6j = -z0*z1*Sq(z2)/(16.0*(z0+z2)*Sq(z0+z2))* x2j_x3j
                    * u1/( Sq(y_u1)*omega_j *u1*t_j1) *u2/( Sq(y_u2)*omega_j *Sq(u2)*Sq(t_j2)) * g_j1 * g_bar_j2 ;
    double term_6k =  z0*z1*Sq(z2)/(16.0*(z1+z2)*Sq(z1+z2))* x2k_x3k
                    * u1/( Sq(y_u1)*omega_k *u1*t_k1) *u2/( Sq(y_u2)*omega_k *Sq(u2)*Sq(t_k2)) * g_k1 * g_bar_k2 ;
    double term_7j = -(z0+z2)*Sq(z2)/(16.0*Sq(z1+z2)) * x2j_x3j 
                    * u1/( Sq(y_u1)*omega_k *u1*t_k1) *u2/( Sq(y_u2)*omega_j *Sq(u2)*Sq(t_j2)) * g_k1 * g_bar_j2 ;
    double term_7k = (z1+z2)*Sq(z2)/(16.0*Sq(z0+z2)) * x2k_x3k
                    * u1/( Sq(y_u1)*omega_j *u1*t_j1) *u2/( Sq(y_u2)*omega_k *Sq(u2)*Sq(t_k2)) * g_j1 * g_bar_k2 ;

    double res = 0.5*Sq(mf) * ( term_1j + term_2j + term_1k + term_2k + term_3a + term_3b + term_4j + term_4k + term_5j + term_5k + term_6j + term_6k + term_7j + term_7k );

    return res;
}




// --------------------------- LONGITUDINAL INTEGRANDS ---------------------------------


double ILdip_massive_LiLogConst(double Q, double z1, double x01sq, double mf) {

    double bessel_inner_fun = sqrt( (Sq(Q)*z1*(1.0-z1) + Sq(mf))*x01sq);
    double dip_res = 0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        dip_res = 0;
    }else{
        dip_res = 4.0*Sq(Q)*Sq(z1)*Sq(1.0-z1)*Sq(gsl_sf_bessel_K0( bessel_inner_fun )) * 
        ( 5.0/2.0 - Sq(M_PI)/3.0 + Sq(log( z1/(1.0-z1) )) + OmegaL_V(Q,z1,mf) + L_dip(Q,z1,mf) );
    }   

    return dip_res;
}

double ILdip_massive_Iab(double Q, double z1, double x01sq, double mf, double xi) {

    double front_factor = 4.0*Sq(Q)*Sq(z1*(1.0-z1));

    double x01 = sqrt( x01sq );

    double kappa_z = sqrt( z1*(1.0-z1)*Sq(Q) + Sq(mf) );
    double bessel_inner_fun = kappa_z * x01;
    double Iab_integrand = 0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        Iab_integrand = 0;
    }else{
        Iab_integrand = gsl_sf_bessel_K0( bessel_inner_fun ) * 1.0/xi * ( -2.0*log(xi)/(1.0-xi) + (1.0+xi)/2.0 ) *
                        (2.0*gsl_sf_bessel_K0( bessel_inner_fun ) - gsl_sf_bessel_K0( sqrt( Sq(kappa_z) + (1.0-z1)*xi/(1.0-xi) *Sq(mf) ) * x01 ) 
                        - gsl_sf_bessel_K0( sqrt( Sq(kappa_z) + z1*xi/(1.0-xi) *Sq(mf) ) * x01 )  );
    }   

    double dip_res = front_factor * Iab_integrand;
    return dip_res;
}

double ILdip_massive_Icd(double Q, double z1, double x01sq, double mf, double xi, double x) {
    // Note that this now uses the earlier parametrization (xi,x) instead of (chi,u). This should be better for numerics.
    double front_factor = 4.0*Sq(Q)*Sq(z1*(1.0-z1));

    double x01 = sqrt( x01sq );

    double kappa_z = sqrt( z1*(1.0-z1)*Sq(Q) + Sq(mf) );
    double bessel_inner_fun = kappa_z * x01;
    double Icd_integrand = 0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        Icd_integrand = 0;
    }else{
        double CLm1 = Sq(z1)*(1.0-xi)/(1.0-z1) * ( -Sq(xi) + x*(1.0-xi)*( 1.0+(1.0-xi)*(1.0+z1*xi/(1.0-z1)) ) / (x*(1.0-xi)+xi/(1.0-z1)) ); // C^L_m(z,x,xi)
        double CLm2 = Sq(1.0-z1)*(1.0-xi)/(z1) * ( -Sq(xi) + x*(1.0-xi)*( 1.0+(1.0-xi)*(1.0+(1.0-z1)*xi/(z1)) ) / (x*(1.0-xi)+xi/(z1)) );; // C^L_m(1-z,x,xi)
        double kappa1 = xi*Sq(mf)/( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(1.0-z1) ) ) * ( xi*(1.0-x) + x*(1.0-z1*(1.0-xi)/(1.0-z1)) ); // kappa(z,x,xi)
        double kappa2 = xi*Sq(mf)/( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(z1) ) ) * ( xi*(1.0-x) + x*(1.0-(1.0-z1)*(1.0-xi)/(z1)) );; // kappa(1-z,x,xi)
        Icd_integrand = gsl_sf_bessel_K0( bessel_inner_fun ) * Sq(mf)*
                     ( (gsl_sf_bessel_K0( bessel_inner_fun ) - gsl_sf_bessel_K0( x01*sqrt( Sq(kappa_z)/(1.0-x) + kappa1 ) ) ) *
                        CLm1 / ( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(1.0-z1) ) * ( x/(1.0-x)*Sq(kappa_z) + kappa1 ) ) +
                        (gsl_sf_bessel_K0( bessel_inner_fun ) - gsl_sf_bessel_K0( x01*sqrt( Sq(kappa_z)/(1.0-x) + kappa2 ) ) ) *
                        CLm2 / ( (1.0-xi)*(1.0-x)*( x*(1.0-xi)+xi/(z1) ) * ( x/(1.0-x)*Sq(kappa_z) + kappa2 ) ));
    }   

    double dip_res = front_factor * Icd_integrand;
    return dip_res;
}

double ILNLOqg_massive_dipole_part_symm(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    // This version takes advantage of the symmetry.

    double front_factor = 4.0*Sq(Q);

    double z0 = 1.0 - z1 - z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double bessel_inner_fun = sqrt( x01sq * ( Sq(Qbar_k) + Sq(mf) ) );


    double facExp = 0.0;
    if(x02sq/x01sq>1e-8 && x02sq/x01sq<5e2){
        facExp = gsl_sf_exp(-x02sq/(gsl_sf_exp(M_EULER)*x01sq));
    }else if (x02sq/x01sq<1e-8){
        facExp = 1.0;
    }else{ facExp = 0.0; }
    

    double res=0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        res = 0;
    }else{
        res = -front_factor * 2.0 * Sq(z1) * ( 2.0*z0*(z0+z2) + Sq(z2) ) * facExp/ x02sq * Sq(gsl_sf_bessel_K0( bessel_inner_fun )) ;
    }
    return res;
}


double ILNLOqg_massive_dipole_part(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    // Here the symmetry argument is not applied.

    double front_factor = 4.0*Sq(Q);

    double z0 = 1.0 - z1 - z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double bessel_inner_fun_k = sqrt( x01sq * ( Sq(Qbar_k) + Sq(mf) ) );
    double bessel_inner_fun_l = sqrt( x01sq * ( Sq(Qbar_l) + Sq(mf) ) );


    double facExp1 = 0.0;
    if(x02sq/x01sq>1e-8 && x02sq/x01sq<5e2){
        facExp1 = gsl_sf_exp(-x02sq/(gsl_sf_exp(M_EULER)*x01sq));
    }else if (x02sq/x01sq<1e-8){
        facExp1 = 1.0;
    }else{ facExp1 = 0.0; }

    double facExp2 = 0.0;
    if(x21sq/x01sq>1e-8 && x21sq/x01sq<5e2){
        facExp2 = gsl_sf_exp(-x21sq/(gsl_sf_exp(M_EULER)*x01sq));
    }else if (x21sq/x01sq<1e-8){
        facExp2 = 1.0;
    }else{ facExp2 = 0.0; }
    

    double res1=0;
    double res2=0;
    if (bessel_inner_fun_k < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        res1 = 0;
    }else{
        res1 = -front_factor * Sq(z1) * ( 2.0*z0*(z0+z2) + Sq(z2) ) * facExp1/ x02sq * Sq(gsl_sf_bessel_K0( bessel_inner_fun_k )) ;
    }


    if (bessel_inner_fun_l < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        res2 = 0;
    }else{
        res2 = -front_factor * Sq(z0) * ( 2.0*z1*(z1+z2) + Sq(z2) ) * facExp2/ x21sq * Sq(gsl_sf_bessel_K0( bessel_inner_fun_l )) ;
    }

    double res = res1 + res2;

    return res;
}



double ILNLOqg_massive_dipole_part_unintegrated(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2) {
    // Version for the qg part that should be integrated over 4 variables from 0 to infinity: t_1, u_1, t_2, u_2.
    // This is the part of qg that is proportional to the dipole amplitude N_01.
    // The variables y_i are to be integrated over 0 to 1. The Jacobian from the change of variables u_i, t_i -> y_ti, y_ui has been taken into account.

    double front_factor = 4.0*Sq(Q);

    double z0 = 1.0 - z1 - z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));

    
    double res1=0;
    double res2=0;

    double jacobian = 1/Sq(y_t1*y_u1*y_t2*y_u2 );
    double t1 = (1.0-y_t1)/y_t1;
    double t2 = (1.0-y_t2)/y_t2;
    double u1 = (1.0-y_u1)/y_u1;
    double u2 = (1.0-y_u2)/y_u2;

    // double intmax = 50.0;
    // double jacobian = intmax*intmax*intmax*intmax;
    // double t1 = intmax*y_t1;
    // double t2 = intmax*y_t2;
    // double u1 = intmax*y_u1;
    // double u2 = intmax*y_u2;


    double gtilde_k_1 = exp( -u1*(Sq(Qbar_k)+Sq(mf)) - x01sq/(4.0*u1) - x02sq / (4.0*t1) );
    double gtilde_k_2 = exp( -u2*(Sq(Qbar_k)+Sq(mf)) - x01sq/(4.0*u2) - x02sq / (4.0*t2) );
    double gtilde_l_1 = exp( -u1*(Sq(Qbar_l)+Sq(mf)) - x01sq/(4.0*u1) - x21sq / (4.0*t1) );
    double gtilde_l_2 = exp( -u2*(Sq(Qbar_l)+Sq(mf)) - x01sq/(4.0*u2) - x21sq / (4.0*t2) );

    res1 = -front_factor * Sq(z1) * ( 2.0*z0*(z0+z2) + Sq(z2) ) * x02sq / 64.0 / (u1*u2*Sq(t1)*Sq(t2)) * exp( -x02sq / ( x01sq * exp( M_EULER ) ) ) * gtilde_k_1 * gtilde_k_2 ;
    
    res2 = -front_factor * Sq(z0) * ( 2.0*z1*(z1+z2) + Sq(z2) ) * x21sq / 64.0 / (u1*u2*Sq(t1)*Sq(t2)) * exp( -x21sq / ( x01sq * exp( M_EULER ) ) ) * gtilde_l_1 * gtilde_l_2 ;

    double res = jacobian * (res1 + res2) ;

    return res;
}



double ILNLOqg_massive_tripole_part_symm(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_u, double y_t) {
    // This version takes advantage of the symmetry.

    double front_factor = 4.0*Sq(Q);

    double z0 = 1.0 - z1 - z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*Sq(z0+z2));
    double omega_l = z1*z2/(z0*Sq(z1+z2));
    double lambda_k = z1*z2/z0;
    double lambda_l = z0*z2/z1;

    double x2_k = sqrt(x02sq);
    double x2_l = sqrt(x21sq);
    double x3_k = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );


    double term1 = 2.0 * Sq(z1) * ( 2.0*z0 * (z0+z2) + Sq(z2) ) * x02sq / 64.0 * 
                    Sq(G_qg( 1, 2, y_u, y_t, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k )) ;
    double term2 = -1.0/16.0 * z1*Sq(z0) * (1.0-z0) * x20x21 * 
                    G_qg( 1, 2, y_u, y_t, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k ) *
                    G_qg( 1, 2, y_u, y_t, Qbar_l, mf, x2_l, x3_l, omega_l, lambda_l );
    double term3 = Sq(mf)/16.0 * Sq(z2)*Sq(z2) * Sq(
                    z1/(z0+z2) * G_qg( 1, 1, y_u, y_t, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k ) -
                    z0/(z1+z2) * G_qg( 1, 1, y_u, y_t, Qbar_l, mf, x2_l, x3_l, omega_l, lambda_l ) );


    double res = front_factor * ( term1 + term2 + term3 );
    return res;
}

double ILNLOqg_massive_tripole_part(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_u, double y_t) {
    // Here the symmetry argument is not applied.

    double front_factor = 4.0*Sq(Q);

    double z0 = 1.0 - z1 - z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*Sq(z0+z2));
    double omega_l = z1*z2/(z0*Sq(z1+z2));
    double lambda_k = z1*z2/z0;
    double lambda_l = z0*z2/z1;

    double x2_k = sqrt(x02sq);
    double x2_l = sqrt(x21sq);
    double x3_k = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );



    double term1k = Sq(z1) * ( 2.0*z0 * (z0+z2) + Sq(z2) ) * x02sq / 64.0 * 
                    Sq(G_qg( 1, 2, y_u, y_t, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k )) ;

    double term1l = Sq(z0) * ( 2.0*z1 * (z1+z2) + Sq(z2) ) * x21sq / 64.0 * 
                    Sq(G_qg( 1, 2, y_u, y_t, Qbar_l, mf, x2_l, x3_l, omega_l, lambda_l )) ;
    double term2 = -1.0/32.0 * z1*z0 * (  z0*(1.0-z0) + z1*(1-z1) )* x20x21 * 
                    G_qg( 1, 2, y_u, y_t, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k ) *
                    G_qg( 1, 2, y_u, y_t, Qbar_l, mf, x2_l, x3_l, omega_l, lambda_l );
    double term3 = Sq(mf)/16.0 * Sq(z2)*Sq(z2) * Sq(
                    z1/(z0+z2) * G_qg( 1, 1, y_u, y_t, Qbar_k, mf, x2_k, x3_k, omega_k, lambda_k ) -
                    z0/(z1+z2) * G_qg( 1, 1, y_u, y_t, Qbar_l, mf, x2_l, x3_l, omega_l, lambda_l ) );


    double res = front_factor * ( term1k + term1l + term2 + term3 );
    return res;
}


double ILNLOqg_massive_tripole_part_unintegrated(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2) {
    // Version for the qg part that should be integrated over 4 variables from 0 to infinity: t_1, u_1, t_2, u_2.
    // This is the part of qg that is proportional to the tripole amplitude N_012.
    // The variables y_i are to be integrated over 0 to 1. The Jacobian from the change of variables u_i, t_i -> y_ti, y_ui has been taken into account.

    double front_factor = 4.0*Sq(Q);

    double z0 = 1.0 - z1 - z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*Sq(z0+z2));
    double omega_l = z1*z2/(z0*Sq(z1+z2));
    double lambda_k = z1*z2/z0;
    double lambda_l = z0*z2/z1;

    double x2_k = sqrt(x02sq);
    double x2_l = sqrt(x21sq);
    double x3_k = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );



    double jacobian = 1/Sq(y_t1*y_u1*y_t2*y_u2 );
    double t1 = (1.0-y_t1)/y_t1;
    double t2 = (1.0-y_t2)/y_t2;
    double u1 = (1.0-y_u1)/y_u1;
    double u2 = (1.0-y_u2)/y_u2;

    // double intmax = 50.0;
    // double jacobian = intmax*intmax*intmax*intmax;
    // double t1 = intmax*y_t1;
    // double t2 = intmax*y_t2;
    // double u1 = intmax*y_u1;
    // double u2 = intmax*y_u2;

    double g_k_1 = exp( -u1*(Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u1) - lambda_k * omega_k * t1* Sq(mf) - x02sq / (4.0*t1) );
    double g_k_2 = exp( -u2*(Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u2) - lambda_k * omega_k * t2* Sq(mf) - x02sq / (4.0*t2) );
    double g_l_1 = exp( -u1*(Sq(Qbar_l)+Sq(mf)) - Sq(x3_l)/(4.0*u1) - lambda_l * omega_l * t1* Sq(mf) - x21sq / (4.0*t1) );
    double g_l_2 = exp( -u2*(Sq(Qbar_l)+Sq(mf)) - Sq(x3_l)/(4.0*u2) - lambda_l * omega_l * t2* Sq(mf) - x21sq / (4.0*t2) );



    double theta_k_1 = (u1 / omega_k - t1 > 0.0) ? 1.0 : 0.0;
    double theta_k_2 = (u2 / omega_k - t2 > 0.0) ? 1.0 : 0.0;
    double theta_l_1 = (u1 / omega_l - t1 > 0.0) ? 1.0 : 0.0;
    double theta_l_2 = (u2 / omega_l - t2 > 0.0) ? 1.0 : 0.0;


    double term1k = theta_k_1 * theta_k_2 *  Sq(z1) * ( 2.0*z0 * (z0+z2) + Sq(z2) ) * x02sq / 64.0 / ( u1*u2*Sq(t1)*Sq(t2) ) * g_k_1 * g_k_2;
    double term1l = theta_l_1 * theta_l_2 * Sq(z0) * ( 2.0*z1 * (z1+z2) + Sq(z2) ) * x21sq / 64.0 / ( u1*u2*Sq(t1)*Sq(t2) ) * g_l_1 * g_l_2;
    double term2 = -theta_k_1 * theta_l_2 * 1.0/32.0 * z1*z0 * (  z0*(1.0-z0) + z1*(1-z1) )* x20x21 / ( u1*u2*Sq(t1)*Sq(t2) ) * g_k_1 * g_l_2 ;
    double term3 = Sq(mf)/16.0 * Sq(z2)*Sq(z2) / (u1*u2*t1*t2) * ( z1/(z0+z2) * theta_k_1 *  g_k_1 - z0/(z1+z2) * theta_l_1 *  g_l_1 ) * ( z1/(z0+z2) * theta_k_2 *  g_k_2 - z0/(z1+z2) * theta_l_2 * g_l_2 ) ;


    double res = front_factor * jacobian *  ( term1k + term1l + term2 + term3 );
    return res;
}



double ILNLOqg_massive_tripole_part_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    // sigma_qg divided into three parts. This part has no additional integrals. Only contains the part proportional to N_012

    double front_factor = 4.0*Sq(Q);

    double z0 = 1-z1-z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*Sq(z0+z2));
    double omega_l = z1*z2/(z0*Sq(z1+z2));

    double x3_k = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );


    double term_k = 1.0/x02sq * Sq(z1) * ( 2.0*z0*(z0+z2) +Sq(z2) ) * Sq( gsl_sf_bessel_K0( sqrt( Sq(Qbar_k) + Sq(mf) ) *sqrt( Sq(x3_k) + omega_k * x02sq ) )  );
    double term_l = 1.0/x21sq * Sq(z0) * ( 2.0*z1*(z1+z2) +Sq(z2) ) * Sq( gsl_sf_bessel_K0( sqrt( Sq(Qbar_l) + Sq(mf) ) *sqrt( Sq(x3_l) + omega_l * x21sq ) )  );
    double term_kl = -2.0 * z0 *z1 * ( z0*(1.0-z0) + z1*(1.0-z1) ) * x20x21 / (x02sq * x21sq) * gsl_sf_bessel_K0( sqrt( Sq(Qbar_k) + Sq(mf)) * sqrt(Sq(x3_k) + omega_k * x02sq) ) * gsl_sf_bessel_K0( sqrt( Sq(Qbar_l) + Sq(mf)) * sqrt(Sq(x3_l) + omega_l * x21sq) );


    double res = front_factor * ( term_k + term_l + term_kl );
    return res;
}



double ILNLOqg_massive_dipole_part_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    // sigma_qg divided into three parts. This part has no additional integrals. Only contains the part proportional to N_01

    double front_factor = 4.0*Sq(Q);
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double z0 = 1-z1-z2;

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*Sq(z0+z2));
    double omega_l = z1*z2/(z0*Sq(z1+z2));
    double lambda_k = z1*z2/z0;
    double lambda_l = z0*z2/z1;

    double x3_k = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );



    double term_k = -1.0/x02sq * Sq(z1) * ( 2.0*z0*(z0+z2) +Sq(z2) ) * exp( -x02sq / x01sq / exp(M_EULER) ) * Sq( gsl_sf_bessel_K0( sqrt( Sq(Qbar_k) + Sq(mf) ) * sqrt(x01sq) )  );
    double term_l = -1.0/x21sq * Sq(z0) * ( 2.0*z1*(z1+z2) +Sq(z2) ) * exp( -x21sq / x01sq / exp(M_EULER) ) * Sq( gsl_sf_bessel_K0( sqrt( Sq(Qbar_l) + Sq(mf) ) * sqrt(x01sq) )  );

    double res = front_factor * ( term_k + term_l );
    return res;
}



double ILNLOqg_massive_tripole_part_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u) {
    // sigma_qg divided into three parts. This part has two additional integrals.

    double front_factor = 4.0*Sq(Q);

    double z0 = 1-z1-z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*Sq(z0+z2));
    double omega_l = z1*z2/(z0*Sq(z1+z2));
    double lambda_k = z1*z2/z0;
    double lambda_l = z0*z2/z1;

    double x3_k = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double u = (1-y_u)/y_u;
    double t_k = y_t * u / omega_k;
    double t_l = y_t * u / omega_l;

    double g_bar_k = exp( -u * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u) ) * exp( -x02sq/(4.0*t_k) ) * ( exp(-t_k*omega_k*lambda_k*Sq(mf)) - 1.0 );
    double g_bar_l = exp( -u * (Sq(Qbar_l)+Sq(mf)) - Sq(x3_l)/(4.0*u) ) * exp( -x21sq/(4.0*t_l) ) * ( exp(-t_l*omega_l*lambda_l*Sq(mf)) - 1.0 );


    double term_k = Sq(z1) * ( 2.0*z0*(z0+z2) +Sq(z2) ) * 1.0/(4.0*Sq(t_k)) * 1/(omega_k*Sq(y_u)) * g_bar_k 
                    * gsl_sf_bessel_K0( sqrt( Sq(Qbar_k) + Sq(mf) ) *sqrt( Sq(x3_k) + omega_k * x02sq )   );
    double term_l = Sq(z0) * ( 2.0*z1*(z1+z2) +Sq(z2) ) * 1.0/(4.0*Sq(t_l)) * 1/(omega_l*Sq(y_u))* g_bar_l 
                    * gsl_sf_bessel_K0( sqrt( Sq(Qbar_l) + Sq(mf) ) *sqrt( Sq(x3_l) + omega_l * x21sq )   );
    double term_kl = -1.0/4.0 * z0*z1 * ( z0*(1.0-z0) + z1*(1.0-z1) ) * x20x21 / Sq(y_u) * (
        1.0/( x21sq * Sq(t_k) * omega_k ) * g_bar_k * gsl_sf_bessel_K0( sqrt( Sq(Qbar_l) + Sq(mf) ) *sqrt( Sq(x3_l) + omega_l * x21sq ))
        + 1.0/( x02sq * Sq(t_l) * omega_l ) * g_bar_l * gsl_sf_bessel_K0( sqrt( Sq(Qbar_k) + Sq(mf) ) *sqrt( Sq(x3_k) + omega_k * x02sq ))
    );

    double res = front_factor * ( term_k + term_l + term_kl );
    return res;
}


double ILNLOqg_massive_tripole_part_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2) {
    // sigma_qg divided into three parts. This part has four additional integrals.

    double front_factor = 4.0*Sq(Q);

    double z0 = 1-z1-z2;
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double Qbar_k = Q*sqrt(z1*(1.0-z1));
    double Qbar_l = Q*sqrt(z0*(1.0-z0));
    double omega_k = z0*z2/(z1*Sq(z0+z2));
    double omega_l = z1*z2/(z0*Sq(z1+z2));
    double lambda_k = z1*z2/z0;
    double lambda_l = z0*z2/z1;

    double x3_k = sqrt( Sq(z0) / Sq(z0+z2) * x02sq + x21sq - 2.0 * z0/(z0+z2) *x20x21 );
    double x3_l = sqrt( Sq(z1) / Sq(z1+z2) * x21sq + x02sq - 2.0 * z1/(z1+z2) *x20x21 );

    double u1 = (1-y_u1)/y_u1;
    double u2 = (1-y_u2)/y_u2;
    double t_k1 = y_t1 * u1 / omega_k;
    double t_l1 = y_t1 * u1 / omega_l;
    double t_k2 = y_t2 * u2 / omega_k;
    double t_l2 = y_t2 * u2 / omega_l;

    double g_bar_k1 = exp( -u1 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u1) ) * exp( -x02sq/(4.0*t_k1) ) * ( exp(-t_k1*omega_k*lambda_k*Sq(mf)) - 1.0 );
    double g_bar_l1 = exp( -u1 * (Sq(Qbar_l)+Sq(mf)) - Sq(x3_l)/(4.0*u1) ) * exp( -x21sq/(4.0*t_l1) ) * ( exp(-t_l1*omega_l*lambda_l*Sq(mf)) - 1.0 );
    double g_bar_k2 = exp( -u2 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u2) ) * exp( -x02sq/(4.0*t_k2) ) * ( exp(-t_k2*omega_k*lambda_k*Sq(mf)) - 1.0 );
    double g_bar_l2 = exp( -u2 * (Sq(Qbar_l)+Sq(mf)) - Sq(x3_l)/(4.0*u2) ) * exp( -x21sq/(4.0*t_l2) ) * ( exp(-t_l2*omega_l*lambda_l*Sq(mf)) - 1.0 );

    double g_k1 = exp( -u1 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u1) ) * exp( -x02sq/(4.0*t_k1) ) * exp(-t_k1*omega_k*lambda_k*Sq(mf));
    double g_l1 = exp( -u1 * (Sq(Qbar_l)+Sq(mf)) - Sq(x3_l)/(4.0*u1) ) * exp( -x21sq/(4.0*t_l1) ) * exp(-t_l1*omega_l*lambda_l*Sq(mf));
    double g_k2 = exp( -u2 * (Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u2) ) * exp( -x02sq/(4.0*t_k2) ) * exp(-t_k2*omega_k*lambda_k*Sq(mf));
    double g_l2 = exp( -u2 * (Sq(Qbar_l)+Sq(mf)) - Sq(x3_l)/(4.0*u2) ) * exp( -x21sq/(4.0*t_l2) ) * exp(-t_l2*omega_l*lambda_l*Sq(mf));


    double term_k = Sq(z1) * ( 2.0*z0*(z0+z2) +Sq(z2) ) * x02sq/64.0 / ( Sq(t_k1*t_k2 * y_u1 * y_u2 * omega_k) ) *g_bar_k1 * g_bar_k2 ;
    double term_l = Sq(z0) * ( 2.0*z1*(z1+z2) +Sq(z2) ) * x21sq/64.0 / ( Sq(t_l1*t_l2 * y_u1 * y_u2 * omega_l) ) *g_bar_l1 * g_bar_l2 ;
    double term_kl = -1.0/32.0 * z1*z0 * ( z1*(1.0-z1) + z0*(1.0-z0) ) * x20x21 / ( Sq( t_k1 * t_l2 * y_u1 * y_u2 ) *omega_k * omega_l ) * g_bar_k1 * g_bar_l2;
    double term_mf = Sq(mf)/16.0 * Sq(z2) * Sq(z2) / ( y_t1 * y_t2 * u1 * u2 * Sq( y_u1 * y_u2 ) ) * (  
        Sq(z1/(z0+z2)) * g_k1 * g_k2 + Sq( z0/(z1+z2) ) * g_l1 * g_l2 - 2.0 * z0/(z1+z2) * z1/(z0+z2) * g_k1 * g_l2
    );


    double res = front_factor * ( term_k + term_l + term_kl + term_mf );
    return res;
}




// --------------------------- TRANSVERSE INTEGRANDS ---------------------------------


double ITdip_massive_0(double Q, double z1, double x01sq, double mf) {

    double x01 = sqrt( x01sq );

    double kappa_z = sqrt( z1*(1.0-z1)*Sq(Q) + Sq(mf) );

    double term1 = Sq( kappa_z * gsl_sf_bessel_K1( x01 *kappa_z ) ) * ( ( Sq(z1) + Sq(1.0-z1) ) * ( 5.0/2.0 - Sq(M_PI)/3.0 + Sq( log(z1/(1.0-z1)) ) + OmegaT_V_unsymmetric(Q, z1, mf) + OmegaT_V_unsymmetric(Q, 1.0-z1, mf) + L_dip(Q,z1,mf)  ) + (2.0*z1-1.0)/2.0 * (OmegaT_N_unsymmetric(Q, z1, mf)-OmegaT_N_unsymmetric(Q, 1.0-z1, mf) ) );
    double term2 = Sq( mf * gsl_sf_bessel_K0( x01 *kappa_z ) ) * ( 3.0 -Sq(M_PI)/3.0 + Sq(log(z1/(1.0-z1))) + OmegaT_V_unsymmetric(Q, z1, mf) + OmegaT_V_unsymmetric(Q, 1.0-z1, mf) + L_dip( Q, z1, mf )  );

    double res= term1 + term2;

    return res;

}


double ITdip_massive_1(double Q, double z1, double x01sq, double mf, double xi)  {
    // One additional integral: xi from 0 to 1.

    double x01 = sqrt( x01sq );

    double kappa_z = sqrt( z1*(1.0-z1)*Sq(Q) + Sq(mf) );

    double term1 = kappa_z * gsl_sf_bessel_K1( x01 * kappa_z) * ( Sq(z1) + Sq(1.0-z1) ) * (IT_V1_unsymmetric(Q, z1, mf, x01, xi) + IT_V1_unsymmetric(Q, 1.0-z1, mf, x01, xi) );
    double term2 = Sq(mf) * gsl_sf_bessel_K0( x01 * kappa_z ) * (IT_VMS1_unsymmetric(Q, z1, mf, x01, xi) + IT_VMS1_unsymmetric(Q, 1.0-z1, mf, x01, xi) );

    double res= term1 + term2;

    return res;

}

double ITdip_massive_2(double Q, double z1, double x01sq, double mf, double y_chi, double y_u) {
    // Two additional integrals: y_chi and y_u, both from 0 to 1


    double x01 = sqrt( x01sq );

    double kappa_z = sqrt( z1*(1.0-z1)*Sq(Q) + Sq(mf) );

    double term1 = kappa_z * gsl_sf_bessel_K1( x01 * kappa_z) * ( 
        ( Sq(z1) + Sq(1.0-z1) ) * (IT_V2_unsymmetric(Q, z1, mf, x01, y_chi, y_u) + IT_V2_unsymmetric(Q, 1.0-z1, mf, x01, y_chi, y_u) )   
        +  (2.0*z1-1.0)/2.0 * (IT_N_unsymmetric(Q, z1, mf, x01, y_chi, y_u) - IT_N_unsymmetric(Q, 1.0-z1, mf, x01, y_chi, y_u) ) );
    double term2 = Sq(mf) * gsl_sf_bessel_K0( x01 * kappa_z ) * (IT_VMS2_unsymmetric(Q, z1, mf, x01, y_chi, y_u) + IT_VMS2_unsymmetric(Q, 1.0-z1, mf, x01, y_chi, y_u) );

    double res= term1 + term2;

    return res;

}


double ITNLOqg_massive_dipole_part_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    return IT_dipole_jk_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq ) + 
           IT_dipole_jkm_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq );
}

double ITNLOqg_massive_tripole_part_I1(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq){
    return IT_tripole_jk_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq ) + 
           IT_tripole_jkm_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq ) + 
           IT_tripole_F_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq ) +
           IT_tripole_Fm_I1( Q, mf, z1, z2, x01sq, x02sq, x21sq );
}

double ITNLOqg_massive_tripole_part_I2(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t, double y_u){
    return IT_tripole_jk_I2( Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t, y_u ) + 
           IT_tripole_jkm_I2( Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t, y_u  ) + 
           IT_tripole_F_I2( Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t, y_u  ) +
           IT_tripole_Fm_I2( Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t, y_u  );
}

double ITNLOqg_massive_tripole_part_I3(double Q, double mf, double z1, double z2, double x01sq, double x02sq, double x21sq, double y_t1, double y_u1, double y_t2, double y_u2) {
    return IT_tripole_jk_I3( Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t1, y_u1, y_t2, y_u2 ) + 
           IT_tripole_jkm_I3( Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t1, y_u1, y_t2, y_u2 ) + 
           IT_tripole_F_I3( Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t1, y_u1, y_t2, y_u2 ) +
           IT_tripole_Fm_I3( Q, mf, z1, z2, x01sq, x02sq, x21sq, y_t1, y_u1, y_t2, y_u2 );
}




