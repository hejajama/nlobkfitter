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


// --------------------------- INTEGRANDS ---------------------------------


double ILdip_massive_LiLogConst(double Q, double z1, double x01sq, double mf) {

    double bessel_inner_fun = sqrt( (Sq(Q)*z1*(1.0-z1) + Sq(mf))*x01sq);
    double dip_res = 0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        dip_res = 0;
    }else{
        dip_res = 4.0*Sq(Q)*Sq(z1)*Sq(1.0-z1)*Sq(gsl_sf_bessel_K0( bessel_inner_fun )) * 
        ( 5.0/2.0 - Sq(M_PI)/6.0 + 1.0/2.0 * Sq(log( z1/(1.0-z1) )) + OmegaL_V(Q,z1,mf) + L_dip(Q,z1,mf) );
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

    double g_k_1 = exp( -u1*(Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u1) - lambda_k * omega_k * t1* Sq(mf) - x02sq / (4.0*t1) );
    double g_k_2 = exp( -u2*(Sq(Qbar_k)+Sq(mf)) - Sq(x3_k)/(4.0*u2) - lambda_k * omega_k * t2* Sq(mf) - x02sq / (4.0*t2) );
    double g_l_1 = exp( -u1*(Sq(Qbar_l)+Sq(mf)) - Sq(x3_l)/(4.0*u1) - lambda_l * omega_l * t1* Sq(mf) - x21sq / (4.0*t1) );
    double g_l_2 = exp( -u2*(Sq(Qbar_l)+Sq(mf)) - Sq(x3_l)/(4.0*u2) - lambda_l * omega_l * t2* Sq(mf) - x21sq / (4.0*t2) );



    double theta_k_1 = u1 / omega_k - t1 ? 1.0 : 0.0;
    double theta_k_2 = u2 / omega_k - t2 ? 1.0 : 0.0;
    double theta_l_1 = u1 / omega_l - t1 ? 1.0 : 0.0;
    double theta_l_2 = u2 / omega_l - t2 ? 1.0 : 0.0;


    double term1k = theta_k_1 * theta_k_2 *  Sq(z1) * ( 2.0*z0 * (z0+z2) + Sq(z2) ) * x02sq / 64.0 / ( u1*u2*Sq(t1)*Sq(t2) ) * g_k_1 * g_k_2;
    double term1l = theta_l_1 * theta_l_2 * Sq(z0) * ( 2.0*z1 * (z1+z2) + Sq(z2) ) * x21sq / 64.0 / ( u1*u2*Sq(t1)*Sq(t2) ) * g_l_1 * g_l_2;
    double term2 = -theta_k_1 * theta_l_2 * 1.0/32.0 * z1*z0 * (  z0*(1.0-z0) + z1*(1-z1) )* x20x21 / ( u1*u2*Sq(t1)*Sq(t2) ) * g_k_1 * g_l_2 ;
    double term3 = Sq(mf)/16.0 * Sq(z2)*Sq(z2) / (u1*u2*t1*t2) * ( z1/(z0+z2) * theta_k_1 *  g_k_1 - z0/(z1+z2) * theta_l_1 *  g_l_1 ) * ( z1/(z0+z2) * theta_k_2 *  g_k_2 - z0/(z1+z2) * theta_l_2 * g_l_2 ) ;


    double res = front_factor * jacobian *  ( term1k + term1l + term2 + term3 );
    return res;
}