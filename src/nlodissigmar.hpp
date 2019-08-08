
#ifndef _NLODIS_SIGMAR_H
#define _NLODIS_SIGMAR_H

#include <string>
#include <vector>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameterState.h>

#include <amplitudelib/amplitudelib.hpp>
#include <tools/interpolation.hpp>
#include "solver.hpp"
#include "ic.hpp"
#include "mv.hpp"
#include "ic_datafile.hpp"
#include "dipole.hpp"
#include "solver.hpp"
#include "nlodis_config.hpp"

#include "cuba-4.2.h"
#include "data.hpp"

#define PARALLEL_CHISQR


using namespace ROOT::Minuit2;
using namespace std;

static inline double Sq(double x){return x*x;}

//namespace sigmar_config{extern double maxy;}

class NLODISFitter : public FCNBase
{
public:
    // MINUIT functions
    double operator() (const vector<double>& par) const; // Calculate Chi^2
    double Up() const {return 1.;} // one-standard-deviation errors for fit parameters.

    // Initialize based on MINUIT parameters
    NLODISFitter(MnUserParameters parameters_);

    // Class configurators
    void AddDataset(Data& d);
    void SetNLO(bool s) { computeNLO = s; }
    void SetSUB(bool s) { UseSub = s; }
    void SetSigma3(bool s) { UseSigma3 = s; }
    void UseImprovedZ2Bound(bool b) { useImprovedZ2Bound = b;}
    void UseConsistentlyBoundLoopTerm(bool b) { useBoundLoop = b;}
    void SetCubaMethod(string s) { cubaMethod = s; }


private:

    bool computeNLO, UseSub, UseSigma3, useImprovedZ2Bound, useBoundLoop;
    string cubaMethod;
    MnUserParameters parameters;
    vector<Data*> datasets;

};


class ComputeSigmaR
{
public:
  typedef  double (ComputeSigmaR::*CmptrMemFn)(double x);
  typedef  double (ComputeSigmaR::*CmptrMemFn_void)(void *x);
  typedef  double (ComputeSigmaR::*z2funpointer)(double x, double q);
  typedef  double (ComputeSigmaR::*subterm_z2upper_fp)(double z);
  typedef  double (ComputeSigmaR::*xrapidity_funpointer)(double x, double q);
  typedef  double (ComputeSigmaR::*xrapidity_NLO_funpointer)(double z2, double z2min, double icX0, double x01sq, double x02sq, double x21sq);
  typedef  double (ComputeSigmaR::*bkkernel_funpointer)(double x01sq, double x02sq, double x21sq);
  typedef  double (ComputeSigmaR::*trbk_rho_funpointer)(double x01sq, double x02sq, double x21sq);
  // old
  typedef  double (ComputeSigmaR::*INLOqg_subterm_fp)(double Q, double x, double z1, double z2,
                                                      double x01sq, double x02sq, double phix0102);


    ComputeSigmaR(AmplitudeLib *ObjectPointer);

    double SigmarLO ( double Q , double xbj, double y) ;
    double SigmarLOmass ( double Q , double xbj, double y, bool charm=false) ;
    double SigmarNLOunsub ( double Q , double xbj, double y) ;
    double SigmarNLOsub ( double Q , double xbj, double y) ;
    double SigmarNLOunsub_UniformZ2Bound ( double Q , double xbj, double y) ;
    double SigmarNLOsub_UniformZ2Bound ( double Q , double xbj, double y) ;
    double SigmarNLOunsub_sigma3 ( double Q , double xbj, double y) ;

    double SigmarNLOunsubRisto ( double Q , double xbj, double y) ;
    double SigmarNLOsubRisto ( double Q , double xbj, double y) ;

    double Structf_LLO ( double Q , double xbj) ;
    double Structf_TLO ( double Q , double xbj) ;
    double Structf_LNLOdip ( double Q , double xbj) ;
    double Structf_TNLOdip ( double Q , double xbj) ;
    double Structf_LNLOdip_z2 ( double Q , double xbj) ;
    double Structf_TNLOdip_z2 ( double Q , double xbj) ;
    double Structf_LNLOqg_sub ( double Q , double xbj) ;
    double Structf_TNLOqg_sub ( double Q , double xbj) ;
    double Structf_LNLOqg_unsub ( double Q , double xbj) ;
    double Structf_TNLOqg_unsub ( double Q , double xbj) ;
    double Structf_LNLOsigma3 ( double Q , double xbj) ;
    double Structf_TNLOsigma3 ( double Q , double xbj) ;

    double Structf_LFULLNLOunsub ( double Q , double xbj) ; // no dipole term in these 4 for now.
    double Structf_TFULLNLOunsub ( double Q , double xbj) ;
    double Structf_LFULLNLOsub ( double Q , double xbj) ;
    double Structf_TFULLNLOsub ( double Q , double xbj) ;

    void SetQuarkMassLight(double qMass_light_){ qMass_light = qMass_light_; }
	void SetQuarkMassCharm(double m) { qMass_charm = m; }
    void SetAlphasScalingC2(double c2_){ alpha_scaling_C2_ = c2_; }
    void SetX0(double x0_){ icX0 = x0_; }
    void SetX0_BK(double x0_){ icX0_bk = x0_; }
    void SetQ0Sqr(double q0_){ icQ0sqr = q0_; }
    void SetY0(double y0_){ icY0 = y0_; }

    void SetRunningCoupling(CmptrMemFn p){AlphabarPTR = p;} // function pointer setter
    void SetRunningCoupling_QG(CmptrMemFn_void p){Alphabar_QG_PTR = p;} // function pointer setter

    void SetImprovedZ2Bound(z2funpointer p){z2limit_PTR = p;} // function pointer setter
    void SetEvolutionX_LO(xrapidity_funpointer p){Xrpdty_LO_PTR = p;}
    void SetEvolutionX_DIP(xrapidity_funpointer p){Xrpdty_DIP_PTR = p;}
    void SetEvolutionX_NLO(xrapidity_NLO_funpointer p){Xrpdty_NLO_PTR = p;}

    void SetSigma3BKKernel(bkkernel_funpointer p){K_kernel_PTR = p;} // function pointer setter
    void SetSubTermKernel(nlodis_config::SubSchemeTermKernel sstk){
        if (sstk == nlodis_config::SUBTERM_LOBK_Z2TOZERO){
            ILNLOqg_subterm_PTR = &ComputeSigmaR::ILNLOqg_subterm_lobk_z2tozero;
            ITNLOqg_subterm_PTR = &ComputeSigmaR::ITNLOqg_subterm_lobk_z2tozero;
        }else if (sstk == nlodis_config::SUBTERM_LOBK_EXPLICIT){
            ILNLOqg_subterm_PTR = &ComputeSigmaR::ILNLOqg_subterm_lobk_explicit;
            ITNLOqg_subterm_PTR = &ComputeSigmaR::ITNLOqg_subterm_lobk_explicit;
        }else if (sstk == nlodis_config::SUBTERM_RESUM){
            ILNLOqg_subterm_PTR = &ComputeSigmaR::ILNLOqg_subterm_resumbk;
            ITNLOqg_subterm_PTR = &ComputeSigmaR::ITNLOqg_subterm_resumbk;
        }else if (sstk == nlodis_config::SUBTERM_KCBK_BEUF){
            ILNLOqg_subterm_PTR = &ComputeSigmaR::ILNLOqg_subterm_kcbk_beuf;
            ITNLOqg_subterm_PTR = &ComputeSigmaR::ITNLOqg_subterm_kcbk_beuf;
        }else if (sstk == nlodis_config::SUBTERM_TRBK_EDMOND){
            ILNLOqg_subterm_PTR = &ComputeSigmaR::ILNLOqg_subterm_trbk_edmond;
            ITNLOqg_subterm_PTR = &ComputeSigmaR::ITNLOqg_subterm_trbk_edmond;
        }
    }
    void SetTRBKRhoPrescription(nlodis_config::TargetRapidityBKRhoPresc rhopresc){
        if (rhopresc == nlodis_config::TRBK_RHO_X_R){rho_PTR = &ComputeSigmaR::rho_rapidity_shift_XR;}
        else if (rhopresc == nlodis_config::TRBK_RHO_MAX_X_Y_R){rho_PTR = &ComputeSigmaR::rho_rapidity_shift_MAX_XYR;}
    }
    void MetaPrescriptionSetter(){
        // NLO: set runningcoupling and C2=Csq for the object.
        ComputeSigmaR::CmptrMemFn alphas_temppointer;
        ComputeSigmaR::CmptrMemFn_void alphas_temppointer_QG;
        if      (nlodis_config::RC_DIS == nlodis_config::DIS_RC_FIXED){
            alphas_temppointer = &ComputeSigmaR::alpha_bar_fixed;
            alphas_temppointer_QG  = &ComputeSigmaR::alpha_bar_QG_fixed;
            cout << "# Using FIXED_LO" << endl;}
        else if (nlodis_config::RC_DIS == nlodis_config::DIS_RC_PARENT){
            alphas_temppointer = &ComputeSigmaR::alpha_bar_running_pd;
            alphas_temppointer_QG  = &ComputeSigmaR::alpha_bar_QG_running_pd;
            cout << "# Using parent dipole RC" << endl;}
        else if (nlodis_config::RC_DIS == nlodis_config::DIS_RC_GUILLAUME){
            alphas_temppointer = &ComputeSigmaR::alpha_bar_running_pd;
            alphas_temppointer_QG  = &ComputeSigmaR::alpha_bar_QG_running_guillaume;
            cout << "# Using Guillaume RC" << endl;}
        else {
            cout << "ERROR: Problem with the choice of runnincoupling. Unkonwn nlodis_config::RC_DIS." << endl;
            exit(1);
            }
        this->SetRunningCoupling(alphas_temppointer);
        this->SetRunningCoupling_QG(alphas_temppointer_QG);

        // Set z2 lower limit settings
        ComputeSigmaR::z2funpointer z2bound_funptr;
        if  (nlodis_config::Z2MINIMUM == nlodis_config::Z2IMPROVED) {z2bound_funptr = &ComputeSigmaR::z2bound_improved;}
        else                             {z2bound_funptr = &ComputeSigmaR::z2bound_simple;}
        this->SetImprovedZ2Bound(z2bound_funptr);

        // Dipole evalution X, (y = log(x0/X))
        // if using 'sub' scheme with z2imp one must use the extended evolution variable for sigma_LO as well.
        ComputeSigmaR::xrapidity_funpointer x_fun_ptr;
        if (nlodis_config::SUB_SCHEME == nlodis_config::SUBTRACTED 
            and nlodis_config::Z2MINIMUM == nlodis_config::Z2IMPROVED) {x_fun_ptr = &ComputeSigmaR::Xrpdty_LO_improved;}
        else {x_fun_ptr = &ComputeSigmaR::Xrpdty_LO_simple;}

        ComputeSigmaR::xrapidity_NLO_funpointer x_nlo_fun_ptr;
        if (config::KINEMATICAL_CONSTRAINT == config::KC_EDMOND_K_MINUS) {x_nlo_fun_ptr = &ComputeSigmaR::Xrpdty_NLO_targetETA;}
        else {x_nlo_fun_ptr = &ComputeSigmaR::Xrpdty_NLO_projectileY;}
        this->SetEvolutionX_LO(x_fun_ptr);
        this->SetEvolutionX_DIP(x_fun_ptr);
        this->SetEvolutionX_NLO(x_nlo_fun_ptr);

        // sigma_3 BK correction
        this->SetSigma3BKKernel(&ComputeSigmaR::K_resum);

        // sub scheme subtraction term choice
        this->SetSubTermKernel(nlodis_config::SUB_TERM_KERNEL);
        
        // trbk rho prescription
        this->SetTRBKRhoPrescription(nlodis_config::TRBK_RHO_PRESC);
    }

    void SetCubaMethod(string s){cubamethod = s;}

	AmplitudeLib* GetDipole() { return ClassScopeDipolePointer; }

//private:
    //variables
    AmplitudeLib *ClassScopeDipolePointer;
    double qMass_light, alpha_scaling_C2_, icX0, icX0_bk, icY0, icQ0sqr;
	double qMass_charm;
    
    // Method and function pointers
    CmptrMemFn AlphabarPTR;
    CmptrMemFn_void Alphabar_QG_PTR;
    z2funpointer z2limit_PTR;
    subterm_z2upper_fp z2upper_qg_PTR;
    xrapidity_funpointer Xrpdty_LO_PTR, Xrpdty_DIP_PTR;
    xrapidity_NLO_funpointer Xrpdty_NLO_PTR;
    bkkernel_funpointer K_kernel_PTR;
    trbk_rho_funpointer rho_PTR;
    INLOqg_subterm_fp ILNLOqg_subterm_PTR, ITNLOqg_subterm_PTR;
    // Cuba integrators
    // old int way: static const int vegas = 1, suave = 2, divonne = 3;
    // int cubamethod;
    string cubamethod;// Currently integrates everything with this. ->? Integrate cumbersome QG AND Bound Dipole Loop contrb. separate one?


    /*
    **  HELPERS & POINTERS
    */
    double Sr(double r, double x);
    double SrY(double r, double y);
    double SrTripole(double x01, double x02, double x21, double x);
    double P(double z);
    double heaviside_theta(double x);

    // running couplings
    double Alphabar(double rsq){ return (this->*AlphabarPTR)(rsq); } // pointer shell function
    double Alphabar_QG(void *userdata){ return (this->*Alphabar_QG_PTR)(userdata); } // pointer shell function
    double alpha_bar_running_pd( double rsq );
    double alpha_bar_running_pd_sharpcutoff( double rsq );
    double alpha_bar_fixed( double rsq );

    // QG term specific couplings
    double alpha_bar_QG_running_pd( void *userdata );
    double alpha_bar_QG_fixed( void *userdata );
    double alpha_bar_QG_running_guillaume( void *userdata );

    // z2 lower bounds
    // pointer, only one since lower bound is set consistently in all terms
    double z2lower_bound( double x, double qsq ) { return (this->*z2limit_PTR)(x,qsq); } // pointer shell function
    double z2bound_simple( double x , double qsq ) { return x/icX0 ; }
    double z2bound_improved( double x , double qsq ) { return (x/icX0)*(icQ0sqr/qsq)  ; }

    // Evolution variables / rapidities
    double Xrpdty_LO( double x, double qsq) { return (this->*Xrpdty_LO_PTR)(x,qsq);}
    double Xrpdty_DIP( double x, double qsq) { return (this->*Xrpdty_DIP_PTR)(x,qsq);}
    double Xrpdty_LO_simple( double x , double qsq ) { return x ; } // X = x0*z2min; z2min = xbj/x0
    double Xrpdty_LO_improved( double x , double qsq ) { return x*icQ0sqr/qsq ; } // X = x0*z2min; z2min = xbj/x0*Q0²/Q²
    double Xrpdty_NLO( double z2, double z2min, double icX0, double x01sq, double x02sq, double x21sq ){
        return (this->*Xrpdty_NLO_PTR)(z2,z2min,icX0,x01sq,x02sq,x21sq);
    };
    double Xrpdty_NLO_projectileY( double z2, double z2min, double icX0, double x01sq, double x02sq, double x21sq );
    double Xrpdty_NLO_targetETA( double z2, double z2min, double icX0, double x01sq, double x02sq, double x21sq );

    // Target rapidity Eta shift calculator(s)
    double rho_rapidity_shift(double x01sq, double x02sq, double x21sq){ return (this->*rho_PTR)(x01sq, x02sq, x21sq) };
    double rho_rapidity_shift_XR(double x01sq, double x02sq, double x21sq);
    double rho_rapidity_shift_MAX_XYR(double x01sq, double x02sq, double x21sq);
    double x_eta_delta_ij_r(double x_ij_sq, double rsq);


    // BK kernels
    double K_kernel( double rsq, double x02sq, double x21sq ) { return (this->*K_kernel_PTR)(rsq,x02sq,x21sq); } // pointer shell function
    double K_resum(double rsq, double x02sq, double x21sq);
    double K_lobk(double rsq, double x02sq, double x21sq);

    // Sub scheme subtraction term pointer and its targets
    double ILNLOqg_subterm(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
        return (this->*ILNLOqg_subterm_PTR)(Q,x,z1,z2,x01sq,x02sq,x21sq);
    }
    double ILNLOqg_subterm_lobk_z2tozero(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ILNLOqg_subterm_lobk_explicit(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ILNLOqg_subterm_resumbk(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ILNLOqg_subterm_kcbk_beuf(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ILNLOqg_subterm_trbk_edmond(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    
    double ITNLOqg_subterm(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
        return (this->*ITNLOqg_subterm_PTR)(Q,x,z1,z2,x01sq,x02sq,x21sq);
    }
    double ITNLOqg_subterm_lobk_z2tozero(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ITNLOqg_subterm_lobk_explicit(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ITNLOqg_subterm_resumbk(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ITNLOqg_subterm_kcbk_beuf(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ITNLOqg_subterm_trbk_edmond(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);

    // Sub scheme subtraction term z2 upper bound pointer and targets
    double z2upper_qg_subterm(double z1){ return (this->*z2upper_qg_PTR)(z1); }   // pointer method
    double z2upper_unity(double z1){ return 1.0; }
    double z2upper_min_z_qqbar(double z1){
        // set z2 integration upper bound to min(z1, 1 - z1)
        double min_momentum;
        if (z1 < 0.5){
            min_momentum = z1;
        } else {
            min_momentum = 1 - z1;
        }
        return min_momentum;
    }


    // cross section integrators
    // Longitudinal
    double ILLO(double Q, double z1, double x01sq) ;
    double ILdip(double Q, double z1, double x01sq) ;
    double Bessel0Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ILNLOqg(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ILNLOqgRisto(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ILNLOsigma3(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;

    double LLOp(double Q, double x) ;
    double LLOpMass(double Q, double x, bool charm) ;
    double LNLOdip(double Q, double x) ;
    double LNLOdip_z2(double Q, double x) ;
    double LNLOqgunsub(double Q, double x) ;
    double LNLOsigma3(double Q, double x) ;
    double LNLOqgsub(double Q, double x) ;
    double LNLOqgunsubRisto(double Q, double x) ;
    double LNLOqgsubRisto(double Q, double x) ;


    // Transverse
    double ITLO(double Q, double z1, double x01sq) ;
    double ITdip(double Q, double z1, double x01sq) ;
    double Bessel1Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ITNLOqg(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ITNLOqgRisto(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ITNLOsigma3(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;


    double TLOp(double Q, double x) ;
    double TLOpMass(double Q, double x, bool charm) ;
    double TNLOdip(double Q, double x) ;
    double TNLOdip_z2(double Q, double x) ;
    double TNLOqgunsub(double Q, double x) ;
    double TNLOsigma3(double Q, double x) ;
    double TNLOqgsub(double Q, double x) ;
    double TNLOqgunsubRisto(double Q, double x) ;
    double TNLOqgsubRisto(double Q, double x) ;
};



void Cuba(string method, int ndim, integrand_t integrand,void *userdata, double *integral, double *error, double *prob);

int integrand_ILLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILdip(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILdip_z2(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgunsub(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILsigma3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgsub(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgunsubRisto(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgsubRisto(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;


int integrand_ITLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITdip(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITdip_z2(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITqgunsub(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITsigma3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITqgsub(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITqgunsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITqgsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;

string PrintVector(vector<double> v);

#endif // _NLODIS_SIGMAR_H
