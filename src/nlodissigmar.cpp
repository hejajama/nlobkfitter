#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
//#include <chrono>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include <amplitudelib/amplitudelib.hpp>
//#include "cuba-4.2.h"
#include "nlodissigmar.hpp"


// HOW TO INTRODUCE THESE?  // TODO
namespace sigmar_config{
  double rmax=25.0; //300.0;
  double rmin=1e-5;
  double Nc=3.0;
  double CF=4.0/3.0;
  double alphaem=1.0/137.0;
  double lambdaqcd=0.241; //GeV
  double Csq=gsl_sf_exp(-2*M_EULER); //TODO
  double sumef=6.0/9.0;
  double X0=1e-2;
  double pifac=1./Sq(2*M_PI);
  double structurefunfac=pifac/alphaem;
}
// REPLACE THIS?            // TODO
using namespace sigmar_config;

///===========================================================================================
//// MAIN FUNCTIONS HERE

// at low Q² Z exchange negligible ( we are only modeling gamma exchange) and the reduced cross section simplifies into
// sigma_r = F_2 - f(y) F_L = F_2 - y^2 / Y_+ * F_L = F_2 - y^2/(1+(1-y)^2) * F_L
double ComputeSigmaR::SigmarLO (double Q , double xbj, double y) {
    //ClassScopeDipolePointer->InitializeInterpolation(xbj);
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LLOp(Q,xbj);
    double FT = fac*TLOp(Q,xbj);
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}
double ComputeSigmaR::SigmarLOmass (double Q , double xbj, double y, bool charm) {
//  ClassScopeDipolePointer->InitializeInterpolation(xbj);
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LLOpMass(Q,xbj,charm);
    double FT = fac*TLOpMass(Q,xbj,charm);
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}


double ComputeSigmaR::SigmarNLO ( double Q , double xbj, double y) {
  //ClassScopeDipolePointer->InitializeInterpolation(xbj);
  double fac = structurefunfac*Sq(Q);
  //cout << "Computing dipole contributions" << endl;
  double FLdip = fac*LNLObeufDIP(Q,xbj);
  double FTdip = fac*TNLObeufDIP(Q,xbj);
  int vegas = 1, suave = 2, divonne = 3;
  int method = vegas;
  //cout << "Computing q->g contributions: L" << endl;
  double FLqg = fac*LNLObeufQGiancu(Q,xbj,method);
  //cout << "Computing q->g contributions: T" << endl;
  double FTqg = fac*TNLObeufQGiancu(Q,xbj,method);
  double FL = FLdip + FLqg;
  double FT = FTdip + FTqg;
  double F2 = FL+FT;
  double fy = Sq(y)/(1+Sq(1-y));
  double sigma = F2 - fy*FL;
  return sigma;
}


///===========================================================================================
double NLODISFitter::operator()(const std::vector<double>& par) const
{
    double chisqr = 0;
    // params to IC
    double qs0sqr       = par[ parameters.Index("qs0sqr")];
    double e_c          = par[ parameters.Index("e_c")];
    double fitsigma0    = par[ parameters.Index("fitsigma0")];
    double alphas_scaling       = par[ parameters.Index("alphascalingC2")]; // MATCH THIS IN IMPACTFACTOR ALPHA_S WHEN NLO
    double anomalous_dimension  = par[ parameters.Index("anomalous_dimension")];
    double initialconditionX0  = par[ parameters.Index("initialconditionX0")];
    double qMass_light  = 0.14; // GeV --- doesn't improve fit at LO
	double qMass_charm = 1.35;
    bool   useMasses    = true;
	bool useCharm = false;


	if (qs0sqr < 0.0001 or qs0sqr > 100 or alphas_scaling < 0.01 or alphas_scaling > 99999 or fitsigma0 < 0.1 or fitsigma0 > 999 or e_c < 1 or e_c > 9999)
		return 9999999;

    cout << "=== Initializing Chi^2 regression === "<< " parameters (" << PrintVector(par) << ")" << endl;

    /*
    // ***Solve resummed BK***
    */
    //cout << "=== Initialize BK solver ===" << endl;

    MV ic;                                            // Initial condition
    ic.SetQsqr(qs0sqr);
    ic.SetAnomalousDimension(anomalous_dimension);
    ic.SetE(e_c);                                     // e_c of MVe parametrization

    Dipole dipole(&ic);
    BKSolver solver(&dipole);
    double maxy = std::log(initialconditionX0/(1e-5)); // divisor=smallest HERA xbj in Q^2 range (1E-05)? //6.8;

    //cout << "=== Solving BK ===" << endl;

    solver.SetAlphasScaling(alphas_scaling);
    solver.Solve(maxy);                                // Solve up to maxy

	dipole.Save("tmp_datafile.dat");

    // Give solution to the AmplitudeLib object
    AmplitudeLib DipoleAmplitude(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());
    DipoleAmplitude.SetX0(initialconditionX0);         // TODO needs to match QG limit X0
	DipoleAmplitude.SetOutOfRangeErrors(false);
    AmplitudeLib *DipolePointer = &DipoleAmplitude;

    ComputeSigmaR SigmaComputer(DipolePointer);
    SigmaComputer.SetX0(initialconditionX0);
    SigmaComputer.SetQuarkMassLight(qMass_light);
	SigmaComputer.SetQuarkMassCharm(qMass_charm);
    // NLO set runnincoupling and C2=Csq for the object.
    //SigmaComputer.SetRunningCoupling(alpha_bar_running_pd); // TODO needs member function pointer? or is there another way?
    SigmaComputer.SetAlphasScalingC2(alphas_scaling);


    cout << "=== Computing Reduced Cross sections ===" << endl;

    /*
     * Loop over datapoints and compute theoretical predictions
     */
    int points=0, totalpoints = 0;
    for (unsigned int dataset=0; dataset<datasets.size(); dataset++)
        totalpoints += datasets[dataset]->NumOfPoints();

    // These loops are trivially parallerizable
    // We only parallerize the inner loop where we have about
    // 250 points (total sigmar) and 50 points (charm)
    for (unsigned int dataset=0; dataset<datasets.size(); dataset++)
    {
#ifdef PARALLEL_CHISQR
//    #pragma omp parallel for schedule(dynamic) reduction(+:chisqr) reduction(+:points)
#endif
        for (int i=0; i<datasets[dataset]->NumOfPoints(); i++)
        {
            double xbj      = datasets[dataset]->xbj(i);
            double y        = datasets[dataset]->y(i);              // inelasticity
            double Q2       = datasets[dataset]->Qsqr(i);
            double Q        = sqrt(Q2);
            double sigmar   = datasets[dataset]->ReducedCrossSection(i);
            double sigmar_err = datasets[dataset]->ReducedCrossSectionError(i);

            double theory=0, theory_charm=0;
            if (!computeNLO && !useMasses) // Compute reduced cross section using leading order impact factors
            {
                theory = (fitsigma0)*SigmaComputer.SigmarLO(Q , xbj , y );
            }
            if (!computeNLO && useMasses)
            {
				theory=0;
				if (datasets[dataset]->OnlyCharm(i)==false)
	                theory = (fitsigma0)*SigmaComputer.SigmarLOmass(Q , xbj , y, false );
				if (xbj*(1.0 + 4.0*1.35*1.35/(Q*Q)) < 0.01 and useCharm)
					theory_charm = (fitsigma0)*SigmaComputer.SigmarLOmass(Q , xbj*(1.0 + 4.0*1.35*1.35/(Q*Q)) , y, true ); // charm
            }

            if (computeNLO) // Full NLO impact factors for reduced cross section
            {
                theory = (fitsigma0)*SigmaComputer.SigmarNLO(Q , xbj , y );
            }

            if (std::isnan(theory) or std::isinf(theory))
            {
                cerr << "Warning: theory result " << theory << " with parameters " << PrintVector(par) << endl;
                theory = 99999999;
            }

            chisqr += datasets[dataset]->Weight()*SQR( (theory+theory_charm - sigmar) / sigmar_err );
            points = points + datasets[dataset]->Weight();

            // Output for plotting
            if(true){
            cout    << setw(10) << xbj          << " "
                    << setw(10) << Q2           << " "
                    << setw(10) << y            << " "
                    << setw(10) << sigmar       << " "
                    << setw(10) << sigmar_err   << " "
                    << setw(10) << theory       << " "
					<< setw(10) << theory_charm << endl;
                  }

        }
    }
    cout << endl << "# Calculated chi^2/N = " << chisqr/points << " (N=" << points << "), parameters (" << PrintVector(par) << ")" << endl<<endl;
    return chisqr;
}

void NLODISFitter::AddDataset(Data& d)
{
    datasets.push_back(&d);
}

NLODISFitter::NLODISFitter(MnUserParameters parameters_)
{
    parameters = parameters_;
    // Muuta INIT? Dipoli solverin alustus tai muuta?
}

string PrintVector(vector<double> v)
{
    stringstream ss;
    for (int i=0; i <v.size(); i++)
    {
        ss << v[i];
        if (i < v.size()-1)
            ss <<", ";
    }
    return ss.str();
}

///===========================================================================================
// CUBA WRAP
namespace cuba_config{
  int method, verbose=0, maxeval=2e7;//1e8;
  double epsrel=1e-2, epsabs=0;
}

void Cuba(int method, int ndim, integrand_t integrand,
	  void *userdata, double *integral, double *error, double *prob) {
  // common arguments
  int ncomp=1, nvec=1, seed=0, mineval=0, last=4;
  int nregions, neval, fail;
  void *spin=NULL;
  char *statefile=NULL;
  if(method==1){
    // Vegas-specific arguments
    int nstart=1000, nincrease=1000, nbatch=1000, gridno=0;
    Vegas(ndim,ncomp,integrand,userdata,nvec,cuba_config::epsrel,
	  cuba_config::epsabs,cuba_config::verbose,seed,mineval,
	  cuba_config::maxeval,nstart,nincrease,nbatch,gridno,statefile,
	  spin,&neval,&fail,integral,error,prob);
  }
  else if(method==2){
    // Suave-specific arguments
    int nnew=10e3, nmin=2;
    double flatness=25; //25;
    Suave(ndim,ncomp,integrand,userdata,nvec,cuba_config::epsrel,
	  cuba_config::epsabs,cuba_config::verbose | last,seed,mineval,
	  cuba_config::maxeval,nnew,nmin,flatness,statefile,spin,
	  &nregions,&neval,&fail,integral,error,prob);
  }
  else if(method==3){
    if(ndim==1) ndim=2;
    // Divonne-specific arguments
    int key1=-4e3, key2=-4e3, key3=1, maxpass=10, ngiven=0, nextra=0;
    double border=1e-4, maxchisq=10, mindeviation=0.1;
    Divonne(ndim,ncomp,integrand,userdata,nvec,cuba_config::epsrel,
	    cuba_config::epsabs,cuba_config::verbose,seed,mineval,
	    cuba_config::maxeval,key1,key2,key3,maxpass,border,maxchisq,
	    mindeviation,ngiven,ndim,NULL,nextra,NULL,statefile,spin,
	    &nregions,&neval,&fail,integral,error,prob);
  }
  else if(method==4){
    if(ndim==1) ndim=2;
    // Cuhre-specific arguments
    int key=0;
    Cuhre(ndim,ncomp,integrand,userdata,nvec,cuba_config::epsrel,
	  cuba_config::epsabs,cuba_config::verbose | last,mineval,
	  cuba_config::maxeval,key,statefile,spin,
	  &nregions,&neval,&fail,integral,error,prob);
  }
}

///===========================================================================================
///===========================================================================================
// CONFIGURATORS

/*void ComputeSigmaR::SetRunningCoupling (double input(double)){
    temp_runningcoupling = input;
}*/

//------------------------------------------------------------------
// HELPERS
double ComputeSigmaR::Sr(double r, double x) {
    double Srx, Nrx;
    if(r<rmin){
        Srx = 1.;
    }else if(r>rmax-1e-7){
        Srx = 0.;
    }else{
        Srx = ClassScopeDipolePointer->S(r,x); //1-Nrx;
    }
    return Srx;
}

double ComputeSigmaR::SrTripole(double x01, double x02, double x21, double x) {
  return Nc/(2*CF)*(Sr(x02,x)*Sr(x21,x) - 1/Sq(Nc)*Sr(x01,x));
}

double ComputeSigmaR::P(double z) {
  return 1.0-z+Sq(z)*0.5;
}

double ComputeSigmaR::heaviside_theta(double x) {
    double result = 0;
    if (x >= 0) result = 1;
    return result;
}

double ComputeSigmaR::temp_runningcoupling( double rsq ) { // TODO POINTER?
  return alpha_bar_running_pd(rsq);
}

double ComputeSigmaR::alpha_bar_running_pd( double rsq ) { // alphabar = Nc/M_PI * alphas
	double logvar = gsl_sf_log(	(4*(alpha_scaling_C2_)/(rsq*Sq(lambdaqcd))	)	);
	double AlphaSres = 12*M_PI/(	(33.0-2*Nc)*logvar );
	if ((logvar<0) || (AlphaSres>0.7)) {
		AlphaSres = 0.7;
		}
    return Nc/M_PI*AlphaSres;
}

///=================================================================
ComputeSigmaR::ComputeSigmaR(AmplitudeLib *ObjectPointer){
    ClassScopeDipolePointer = ObjectPointer;
}

struct Userdata{
    ComputeSigmaR* ComputerPtr;
    double Q, xbj;
    double qMass_light;
    double icX0;
};


///===========================================================================================
// INTEGRATORS for cross sections

/*
// --- L L L --- LO -------- L L L --- LO -------- L L L --- LO -----------
*/
double ComputeSigmaR::ILLO(double Q, double z1, double x01sq) {
    double res;
    if(Q*sqrt(z1*(1.0-z1)*x01sq)>1e-7&&Q*sqrt(z1*(1.0-z1)*x01sq)<1e2){
        res = 4.0*Sq(Q)*Sq(z1)*Sq(1.0-z1)*Sq(gsl_sf_bessel_K0(Q*sqrt(z1*(1.0-z1)*x01sq)));
    }else{
        res = 0;
    }
    return res;
}

// MASSLESS
int integrand_ILLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=rmax*x[1];
    double x01sq=x01*x01;

    double res=(1.0-(Optr->Sr(x01,xbj)))*(Optr->ILLO(Q,z1,x01sq))*x01;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::LLOp(double Q, double x ) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.ComputerPtr=this;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Cuba(2,ndim,integrand_ILLOp,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*rmax*integral;
}

// MASS
int integrand_ILLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double qmass=dataptr->qMass_light;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=rmax*x[1];
    //double x01sq=x01*x01;

    double af = sqrt( Sq(Q)*z1*(1-z1) + Sq(qmass) );
    double impactfac = 4.0*Sq(Q*(z1)*(1.0-z1)*gsl_sf_bessel_K0(af*x01));
    double res=(1.0-(Optr->Sr(x01,xbj)))*(impactfac)*x01;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::LLOpMass(double Q, double x, bool charm) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
	if (charm==false)
    	userdata.qMass_light=qMass_light;
	else
		userdata.qMass_light=qMass_charm;
    userdata.ComputerPtr=this;
	double ef = sumef;
	if (charm) ef = 4.0/9.0;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*ef;
    Cuba(2,ndim,integrand_ILLOpMass,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*rmax*integral;
}


/*
// --- T T T --- LO -------- T T T --- LO -------- T T T --- LO -----------
*/

double ComputeSigmaR::ITLO(double Q, double z1, double x01sq) {
  double res;
  if(Q*sqrt(z1*(1.0-z1)*x01sq)>1e-7&&Q*sqrt(z1*(1.0-z1)*x01sq)<1e2){
    res = Sq(Q)*z1*(1.0-z1)*(1.0-2.0*z1+2.0*Sq(z1))*Sq(gsl_sf_bessel_K1(Q*sqrt(z1*(1.0-z1)*x01sq)));
  }else{
    res = 0;
  }
  return res;
}

// MASSLESS
int integrand_ITLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
  Userdata *dataptr = (Userdata*)userdata;
  double Q=dataptr->Q;
  double xbj=dataptr->xbj;
  ComputeSigmaR *Optr = dataptr->ComputerPtr;
  double z1=x[0];
  double x01=rmax*x[1];
  double x01sq=x01*x01;

  double res=(1.0-(Optr->Sr(x01,xbj)))*(Optr->ITLO(Q,z1,x01sq))*x01;
  if(gsl_finite(res)==1){
    *f=res;
  }else{
    *f=0;
  }
  return 0;
}

double ComputeSigmaR::TLOp(double Q, double x) {
  double integral, error, prob;
  const int ndim=2;
  Userdata userdata;
  userdata.Q=Q;
  userdata.xbj=x;
  userdata.ComputerPtr=this;
  double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
  Cuba(2,ndim,integrand_ITLOp,&userdata,&integral,&error,&prob);
  return fac*2.0*M_PI*rmax*integral;
}

// MASS
int integrand_ITLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
  Userdata *dataptr = (Userdata*)userdata;
  double Q=dataptr->Q;
  double xbj=dataptr->xbj;
  double qmass=dataptr->qMass_light;
  ComputeSigmaR *Optr = dataptr->ComputerPtr;
  double z1=x[0];
  double x01=rmax*x[1];
  //double x01sq=x01*x01;

  double af = sqrt( Sq(Q)*z1*(1.0-z1) + Sq(qmass) );
  double impactfac = (1.0-2.0*z1+2.0*Sq(z1))*Sq(af*gsl_sf_bessel_K1(af*x01)) + Sq( qmass*gsl_sf_bessel_K0( af*x01 ) );
  double res=(1.0-(Optr->Sr(x01,xbj)))*(impactfac)*x01;
  if(gsl_finite(res)==1){
    *f=res;
  }else{
    *f=0;
  }
  return 0;
}

double ComputeSigmaR::TLOpMass(double Q, double x, bool charm) {
  double integral, error, prob;
  const int ndim=2;
  Userdata userdata;
  userdata.Q=Q;
  userdata.xbj=x;
  if (charm==false)
  	userdata.qMass_light=qMass_light;
  else
    userdata.qMass_light = qMass_charm;
  userdata.ComputerPtr=this;

  double ef = sumef;
  if (charm) ef = 4.0/9.0;

  double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
  Cuba(2,ndim,integrand_ITLOpMass,&userdata,&integral,&error,&prob);
  return fac*2.0*M_PI*rmax*integral;
}



/*
// --- L L L --- NLO -------- L L L --- NLO -------- L L L --- NLO -----------
*/

///*
// NLO
double ComputeSigmaR::ILbeufDIP(double Q, double z1, double x01sq) {
    double fac1 = Sq(z1)*Sq(1.0 - z1);
    double facNLO;

    if(Q*sqrt(x01sq)>1e-7&&Q*sqrt(x01sq)<1e3){
        facNLO = 4.0*Sq(Q*gsl_sf_bessel_K0(Q*sqrt(z1*(1-z1)*x01sq)));
    }else{
        facNLO = 0;
        }
    double res = fac1*facNLO;
    return res;
}

int integrand_ILbeufDIP(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=rmax*x[1];
    double x01sq=Sq(x01);

    double alphabar=Optr->temp_runningcoupling(x01sq); //2;
    double alphfac=alphabar*CF/Nc;
    double SKernel = 1.0 - Optr->Sr(x01,xbj);
    double regconst = 5.0/2.0 - Sq(M_PI)/6.0;
    double res;

    res = SKernel*(Optr->ILbeufDIP(Q,z1,x01sq))*x01*( 1.0 + alphfac*(0.5*Sq(log(z1/(1-z1))) + regconst ) );
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::Bessel0Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq, double X3sq) {
		//double X3sq = z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq;
		double facNLO;
		double x01=sqrt(x01sq);
		double x02=sqrt(x02sq);
		double x21=sqrt(x21sq);

		if(Q*sqrt(X3sq)>1e-7&&Q*sqrt(X3sq)<1e2){
			facNLO = 4.0*Sq(Q*gsl_sf_bessel_K0(Q*sqrt(X3sq)));
		}else{
			facNLO = 0;
		}
		double res = facNLO*(1-SrTripole(x01,x02,x21,x));
		return res;
}

double ComputeSigmaR::ILNLObeufQG(double Q, double x, double z1, double z2, double x01sq, double x02sq, double phix0102) {
  	double x21sq = x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
  	double x20x21 = -0.5*(x01sq - x21sq - x02sq);
    double X3sq 	= z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq;
    double X010sq = z1*(1.0 - z1)*x01sq;

  	double fac1 = Sq(z1)*Sq(1.0 - z1);
    double xi 	= z2/(1.0-z1);
  	double fun1 = 1+Sq(1-xi);
    double fun2 = x20x21/(x02sq*x21sq);
  	double fac2 = 1/x02sq - fun2;
    double fac3 = Sq(xi)*fun2;

    double facNLO1 = Bessel0Tripole(Q, x, z1, z2, x01sq, x02sq, x21sq, X3sq);
    double facNLO2 = Bessel0Tripole(Q, x, z1, z2, x01sq, 0    , x01sq, X010sq);

  	double res = fac1*((fun1*fac2)*(facNLO1 - facNLO2) + fac3*facNLO1 );
  	return res;
}

int integrand_ILbeufQGiancu(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double x0lim=xbj/X0;
    double z1=(1.0-x0lim)*x[0];
    double z2=(1.0-z1-x0lim)*x[1]+x0lim;
    double x01=rmax*x[2];
    double x02=rmax*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double jac=(1.0-x0lim)*(1.0-z1-x0lim);
    double Xrpdt=xbj/z2;

    double alphabar=Optr->temp_runningcoupling(x01sq);
    double alphfac=alphabar*CF/Nc;

    double res =   jac*alphfac*( Optr->ILNLObeufQG(Q,Xrpdt,z1,z2,x01sq,x02sq,phix0102) )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}


// integrations
double ComputeSigmaR::LNLObeufDIP(double Q, double x) {
    double integral, error, prob;
    const int ndim=2;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.ComputerPtr=this;
    Cuba(2,ndim,integrand_ILbeufDIP,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*rmax*integral;
}

double ComputeSigmaR::LNLObeufQGiancu(double Q, double x, int mthd) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(mthd,ndim,integrand_ILbeufQGiancu,&userdata,&integral,&error,&prob);
    return 2*fac*2.0*M_PI*rmax*rmax*integral;
}
//*/



/*
// --- T T T --- NLO -------- T T T --- NLO -------- T T T --- NLO -----------
*/

///*
// NLO
double ComputeSigmaR::ITbeufDIP(double Q, double z1, double x01sq) {
    double fac1 = z1*(1-z1)*(Sq(z1)+Sq(1.0 - z1));
    double facNLO;

    if(Q*sqrt(x01sq)>1e-7&&Q*sqrt(x01sq)<1e3){
        facNLO = Sq(Q*gsl_sf_bessel_K1(Q*sqrt(z1*(1-z1)*x01sq)));
    }else{
        facNLO = 0;
        }
    double res = fac1*facNLO;
    return res;
}

int integrand_ITbeufDIP(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=rmax*x[1];
    double x01sq=Sq(x01);

    double alphabar=Optr->temp_runningcoupling(x01sq);
    double alphfac=alphabar*CF/Nc;
    double SKernel = 1.0 - Optr->Sr(x01,xbj);
    double regconst = 5.0/2.0 - Sq(M_PI)/6.0;
    double res;

    res = SKernel*( Optr->ITbeufDIP(Q,z1,x01sq) )*x01*( 1.0 + alphfac*(0.5*Sq(log(z1/(1-z1))) + regconst ) );
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::Bessel1Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq, double X3sq) {
		//double X3sq = z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq;
		double facNLO;
		double x01=sqrt(x01sq);
		double x02=sqrt(x02sq);
		double x21=sqrt(x21sq);

		if(Q*sqrt(X3sq)>1e-7&&Q*sqrt(X3sq)<1e2){
			facNLO = Sq(Q*gsl_sf_bessel_K1(Q*sqrt(X3sq)));
		}else{
			facNLO = 0;
		}
		double res = facNLO*(1-SrTripole(x01,x02,x21,x));
		return res;
}

double ComputeSigmaR::ITNLObeufQG(double Q, double x, double z1, double z2, double x01sq, double x02sq, double phix0102) {
    double x21sq 	= x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);
    double X3sq 	= z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq;
    double X010sq = z1*(1.0 - z1)*x01sq;

    double fac1 	= z1*(1.0 - z1);
    double xi 		= z2/(1.0-z1);
    double fac11	= Sq(z1)+Sq(1-z1);
    double fac12 	= 1+Sq(1-xi);
    double fun2 	= x20x21/(x02sq*x21sq);
    double fac13 	= 1/x02sq - fun2;
    double fac3 	= Sq(xi)*(	fac11*fun2	+ 2*fac1*(1-xi)*x20x21/(x02sq*X3sq)	-	(1-z1)*(1-xi)*(z1+xi-z1*xi)/X3sq	);

    double facNLO1 = Bessel1Tripole(Q, x, z1, z2, x01sq, x02sq, x21sq, X3sq);
    double facNLO2 = Bessel1Tripole(Q, x, z1, z2, x01sq, 0    , x01sq, X010sq);

    double res = fac1*((fac11*fac12*fac13)*(facNLO1 - facNLO2) + fac3*facNLO1);
    return res;
}

int integrand_ITbeufQGian(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double x0lim=xbj/X0;
    double z1=(1.0-x0lim)*x[0];
    double z2=(1.0-z1-x0lim)*x[1]+x0lim;
    double x01=rmax*x[2];
    double x02=rmax*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double jac=(1.0-x0lim)*(1.0-z1-x0lim);
    double Xrpdt=xbj/z2;

    double alphabar=Optr->temp_runningcoupling(x01sq);
    double alphfac=alphabar*CF/Nc;

    double res = jac*alphfac*( Optr->ITNLObeufQG(Q,Xrpdt,z1,z2,x01sq,x02sq,phix0102) )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

// integrations
double ComputeSigmaR::TNLObeufDIP(double Q, double x) {
    double integral, error, prob;
    const int ndim=2;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.ComputerPtr=this;
    Cuba(2,ndim,integrand_ITbeufDIP,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*rmax*integral;
}

double ComputeSigmaR::TNLObeufQGiancu(double Q, double x, int mthd) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(mthd,ndim,integrand_ITbeufQGian,&userdata,&integral,&error,&prob);
    return 2*fac*2.0*M_PI*rmax*rmax*integral;
}
//*/
