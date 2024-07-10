#ifndef _INVDIP_H
#define _INVDIP_H

#include <string>
#include <vector>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameterState.h>

#include "cuba-4.2.h"
#include "data.hpp"

#define PARALLEL_CHISQR

using namespace ROOT::Minuit2;
using namespace std;

static inline double Sq(double x){return x*x;}

class InverseDipoleFitter : public FCNBase
{
public:
    // MINUIT functions
    double operator() (const vector<double>& par) const; // Calculate Chi^2
    double Up() const {return 1.;} // one-standard-deviation errors for fit parameters.

    // Initialize based on MINUIT parameters
    InverseDipoleFitter(MnUserParameters parameters_);

    // Class configurators
    void AddDataset(Data& d);
    void AddDipoleGrid(std::vector< std::tuple<double, double, double> > v){private_dipoleGrid = v;}
    void SetNLO(bool s) { computeNLO = s; }
    void SetSUB(bool s) { UseSub = s; }
    void SetSigma3(bool s) { UseSigma3 = s; }
    void UseImprovedZ2Bound(bool b) { useImprovedZ2Bound = b;}
    void UseConsistentlyBoundLoopTerm(bool b) { useBoundLoop = b;}
    void SetCubaMethod(string s) { cubaMethod = s; }


private:

    bool computeNLO, UseSub, UseSigma3, useImprovedZ2Bound, useBoundLoop;
    std::vector< std::tuple<double, double, double> > private_dipoleGrid;
    string cubaMethod;
    MnUserParameters parameters;
    vector<Data*> datasets;

};
// Evolution variables / rapidities
double Xrpdty_LO( double x , double qsq ) { return x ; } // X = x0*z2min; z2min = xbj/x0 // This is used unless SUB scheme is used at NLO

double Structf_LLO ( double Q , double xbj) ;
double Structf_TLO ( double Q , double xbj) ;
double LLOp(double Q, double x) ;
double TLOp(double Q, double x) ;

void Cuba(string method, int ndim, integrand_t integrand,void *userdata, double *integral, double *error, double *prob);
string PrintVector(vector<double> v);

int integrand_ILLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;

int integrand_ITLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;

#endif
