
#ifndef _NLODIS_SIGMAR_H
#define _NLODIS_SIGMAR_H

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

//#define PARALLEL_CHISQR

using namespace ROOT::Minuit2;
using namespace std;

static inline double Sq(double x){return x*x;}

//namespace sigmar_config{extern double maxy;}

class NLODISFitter : public FCNBase
{
public:
    // MINUIT functions
    double operator() (const vector<double>& par) const; // Calculate Chi^2
    double Up() const {return 1.;} // MINUIT default possibly 1.0?

    // Initialize based on MINUIT parameters
    NLODISFitter(MnUserParameters parameters_);

    // Class configurators
    void AddDataset(Data& d);
    void SetNLO(bool s) { computeNLO = s; }

private:

    bool computeNLO;
    MnUserParameters parameters;
    vector<Data*> datasets;

};

class ComputeSigmaR
{
public:

    ComputeSigmaR(AmplitudeLib *ObjectPointer);

    double SigmarLO ( double Q , double xbj, double y) ;
    double SigmarLOmass ( double Q , double xbj, double y, bool charm=false) ;
    double SigmarNLO ( double Q , double xbj, double y) ;

    void SetQuarkMassLight(double qMass_light_){ qMass_light = qMass_light_; }
	void SetQuarkMassCharm(double m) { qMass_charm = m; }
    void SetAlphasScalingC2(double c2_){ alpha_scaling_C2_ = c2_; }
    void SetX0(double x0_){ icX0 = x0_; }
    //void SetRunningCoupling (double input(double));

	AmplitudeLib* GetDipole() { return ClassScopeDipolePointer; }

//private:
    //variables
    AmplitudeLib *ClassScopeDipolePointer;
    //double (*temp_runningcoupling)(double); // TODO Pointer to a running coupling function.
    double temp_runningcoupling(double);
    double qMass_light, alpha_scaling_C2_, icX0;
	double qMass_charm;

    // helpers
    double Sr(double r, double x);
    double SrTripole(double x01, double x02, double x21, double x);
    double P(double z);
    double heaviside_theta(double x);
    double alpha_bar_running_pd( double rsq );

    // cross section integrators
    // Longitudinal
    double ILLO(double Q, double z1, double x01sq) ;
    double ILbeufDIP(double Q, double z1, double x01sq) ;
    double Bessel0Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq, double X3sq) ;
    double ILNLObeufQG(double Q, double x, double z1, double z2, double x01sq, double x02sq, double phix0102) ;

    double LLOp(double Q, double x) ;
    double LNLObeufDIP(double Q, double x) ;
    double LNLObeufQGiancu(double Q, double x, int mthd) ;
    double LLOpMass(double Q, double x, bool charm=false) ;


    // Transverse
    double ITLO(double Q, double z1, double x01sq) ;
    double ITbeufDIP(double Q, double z1, double x01sq) ;
    double Bessel1Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq, double X3sq) ;
    double ITNLObeufQG(double Q, double x, double z1, double z2, double x01sq, double x02sq, double phix0102) ;

    double TLOp(double Q, double x) ;
    double TNLObeufDIP(double Q, double x) ;
    double TNLObeufQGiancu(double Q, double x, int mthd) ;
    double TLOpMass(double Q, double x, bool charm=false) ;

};



void Cuba(int method, int ndim, integrand_t integrand,void *userdata, double *integral, double *error, double *prob);

int integrand_ILLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILbeufDIP(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILbeufQGiancu(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;


int integrand_ITLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITbeufDIP(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITbeufQGian(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;

string PrintVector(vector<double> v);

#endif // _NLODIS_SIGMAR_H
