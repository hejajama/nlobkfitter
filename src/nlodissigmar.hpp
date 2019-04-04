
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
    double Up() const {return 2.;} // MINUIT default possibly 1.0?

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
  typedef  double (ComputeSigmaR::*bkkernel_funpointer)(double x01sq, double x02sq, double x21sq);

    ComputeSigmaR(AmplitudeLib *ObjectPointer);

    double SigmarLO ( double Q , double xbj, double y) ;
    double SigmarLOmass ( double Q , double xbj, double y) ;
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
    void SetAlphasScalingC2(double c2_){ alpha_scaling_C2_ = c2_; }
    void SetX0(double x0_){ icX0 = x0_; }
    void SetQ0Sqr(double q0_){ icQ0sqr = q0_; }
    void SetY0(double y0_){ icY0 = y0_; }

    void SetRunningCoupling(CmptrMemFn p){AlphabarPTR = p;} // function pointer setter
    void SetRunningCoupling_QG(CmptrMemFn_void p){Alphabar_QG_PTR = p;} // function pointer setter

    void SetImprovedZ2Bound(z2funpointer p){z2limit_PTR = p;} // function pointer setter
    void SetSigma3BKKernel(bkkernel_funpointer p){K_kernel_PTR = p;} // function pointer setter

    void SetCubaMethod(string s){cubamethod = s;}


//private:
    //variables
    AmplitudeLib *ClassScopeDipolePointer;
    double qMass_light, alpha_scaling_C2_, icX0, icY0, icQ0sqr;
    CmptrMemFn AlphabarPTR;
    CmptrMemFn_void Alphabar_QG_PTR;
    z2funpointer z2limit_PTR;
    bkkernel_funpointer K_kernel_PTR;
    // Cuba integrators
    // old int way: static const int vegas = 1, suave = 2, divonne = 3;
    // int cubamethod;
    string cubamethod;// Currently integrates everything with this. ->? Integrate cumbersome QG AND Bound Dipole Loop contrb. separate one?

    // helpers
    double Sr(double r, double x);
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
    double z2lower_bound( double x, double qsq ) {return (this->*z2limit_PTR)(x,qsq); } // pointer shell function
    double z2bound_simple( double x , double qsq ) { return x/icX0 ; }
    double z2bound_improved( double x , double qsq ) { return (x/icX0)*(icQ0sqr/qsq)  ; }
    // BK kernels
    double K_kernel( double rsq, double x02sq, double x21sq ) {return (this->*K_kernel_PTR)(rsq,x02sq,x21sq); } // pointer shell function
    double K_resum (double rsq, double x02sq, double x21sq);
    double K_lobk (double rsq, double x02sq, double x21sq);




    // cross section integrators
    // Longitudinal
    double ILLO(double Q, double z1, double x01sq) ;
    double ILdip(double Q, double z1, double x01sq) ;
    double Bessel0Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq, double X3sq) ;
    double ILNLOqg(double Q, double x, double z1, double z2, double x01sq, double x02sq, double phix0102) ;
    double ILNLOqgRisto(double Q, double x, double z1, double z2, double x01sq, double x02sq, double phix0102) ;
    double ILNLOsigma3(double Q, double x, double z1, double z2, double x01sq, double x02sq, double phix0102) ;

    double LLOp(double Q, double x) ;
    double LLOpMass(double Q, double x) ;
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
    double Bessel1Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq, double X3sq) ;
    double ITNLOqg(double Q, double x, double z1, double z2, double x01sq, double x02sq, double phix0102) ;
    double ITNLOqgRisto(double Q, double x, double z1, double z2, double x01sq, double x02sq, double phix0102) ;
    double ITNLOsigma3(double Q, double x, double z1, double z2, double x01sq, double x02sq, double phix0102) ;


    double TLOp(double Q, double x) ;
    double TLOpMass(double Q, double x) ;
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
