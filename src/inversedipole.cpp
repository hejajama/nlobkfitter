#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <csignal>
#include <ctime>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <tools/tools.hpp>
#include <sstream>
#include <unistd.h>

#include <amplitudelib/amplitudelib.hpp>
#include <tools/interpolation.hpp>
#include "solver.hpp"
#include "ic.hpp"
#include "mv.hpp"
#include "ic_datafile.hpp"
#include "dipole.hpp"
#include "solver.hpp"

#include "data.hpp"
// #include "nlodissigmar.hpp"
#include "inversedipole.hpp"
#include "nlodis_config.hpp"
#include "helper.hpp"

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnScan.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>

 #include <gsl/gsl_errno.h>

using namespace std;
using namespace ROOT::Minuit2;

namespace invdip_config{
    // double rmax=30.0; //300.0;
    // double rmin=1e-5;
    double Nc=3.0;
    double sumef=6.0/9.0; // light quarks uds only.
    double CF=4.0/3.0; // (Nc()*Nc()-1.0)/(2.0*Nc());
    double alphaem=1.0/137.0;
    double lambdaqcd=0.241; //GeV
    double structurefunfac=1./(Sq(2*M_PI)*alphaem);
    double icX0_bk=1e-2;
    AmplitudeLib* DipolePointer;
    std::vector< std::tuple<double, double, double> > dipoleGrid;
    std::vector< double > rvals, xvals, Svals;
    string cubaMethod = "vegas";
}
using namespace invdip_config;



//void ErrHandler(const char * reason,const char * file,int line,int gsl_errno);
int errors_mmyiss2;
void ErrHandlerCustom2(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{   // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    // 15: underflows

    if (gsl_errno == 11) return; // ignore max subdivision errors
    if (gsl_errno == 18) return; // roundoff errors from small r?
    if (gsl_errno == 15) return; // underflow // safe to ignore?
    // if (gsl_errno == 16) return; // overflow
    // Ugly hack, comes from the edges of the z integral in virtual_photon.cpp
    // Overflows come from IPsat::bint when it is done analytically
    // Hope is that these errors are handled correctly everywhere
    errors_mmyiss2++;
    std::cerr << file << ":"<< line <<": Error " << errors_mmyiss2 << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}

double Sr_read_grid(double r, double x){
    for (auto elem : dipoleGrid){
        auto [ri, xi, Si] = elem;
        if (r == ri && x == xi) return Si;
    }
    cout << "r and x values not in grid!" << r << ", " << x << ";" << endl;
    exit(0);
    return 0;
}

double Sr_interpolator(double r, double x){
    double Srx_interp;
    // find nearest element values
    double S_bb, S_ab, S_ba, S_aa;
    // auto r_above = std::upper_bound(rvals.begin(), rvals.end(), r);
    // auto x_above = std::upper_bound(xvals.begin(), xvals.end(), r);
    // int ri_above = std::distance(rvals.begin(), r_above);
    // int ri_below = ri_above-1;
    // int xi_above = std::distance(xvals.begin(), x_above);
    // int xi_below = xi_above-1;

    double r1_;
    double r2_;
    double y1_;
    double y2_;
    // cout << endl << "asking r " << r << " x " << x << endl;
    // S_bb = Sr_read_grid(r1_, y1_);
    // S_ab = Sr_read_grid(r2_, y1_);
    // S_ba = Sr_read_grid(r1_, y2_);
    // S_aa = Sr_read_grid(r2_, y2_);

    for(int i=0; i<rvals.size(); i++){
        if ((rvals[i] < r && r <= rvals[i+1])){
            r1_ = rvals[i];
            r2_ = rvals[i+1];
        }
    }
    for(int i=0; i<xvals.size(); i++){
        if ((xvals[i] >= x && x > xvals[i+1])){
            y1_ = xvals[i+1];
            y2_ = xvals[i];
        }
    }
    // cout << "r1 " << r1_ << " r2 " << r2_ << "y1 " << y1_ << " y2 " << y2_  << endl;
    S_bb = Sr_read_grid(r1_, y1_);
    S_ab = Sr_read_grid(r2_, y1_);
    S_ba = Sr_read_grid(r1_, y2_);
    S_aa = Sr_read_grid(r2_, y2_);

    // LINEAR_LINEAR
    // int rind2 = rind+1;
    // int yind2 = yind+1;
    // if (rind2 >= rvals.size()) return 1.0;
    // if (yind2 >= yvals.size()) return n[yind][rind]; // this should never happen

    double y = x;
    // double bilin =  (1.0/( (r2_ - r1_)*(y2_ - y1_) ))*( (n[yind][rind])*(r2_ - r)*(y2_ - y) + (n[yind][rind2])*(r - r1_)*(y2_ - y) + (n[yind2][rind])*(r2_ - r)*(y - y1_) + (n[yind2][rind2])*(r - r1_)*(y - y1_) );
    double bilin =  (1.0/( (r2_ - r1_)*(y2_ - y1_) ))*( (S_bb)*(r2_ - r)*(y2_ - y) + (S_ab)*(r - r1_)*(y2_ - y) + (S_ba)*(r2_ - r)*(y - y1_) + (S_aa)*(r - r1_)*(y - y1_) );
    return bilin;

    return Srx_interp;
}

double Sr(double r, double x){
    double Srx;//, Nrx;
    // cout << "Sr r x: " << r << " " << x << endl;
    if(r<nlodis_config::MINR){
        Srx = 1.;
    }else if(r>nlodis_config::MAXR-1e-7){
        Srx = 0.;
    }else{

        if (x > icX0_bk){
            Srx = Sr_interpolator(r,icX0_bk) ; //1-Nrx;
        } else if (x < xvals.back()){
            Srx = Sr_interpolator(r,xvals.back()) ; //1-Nrx;
        } else {
            Srx = Sr_interpolator(r,x) ; //1-Nrx;
        }
    }
    return Srx;
}

/*
// --- Structure functions -----------
*/
double SigmarLO (double Q , double xbj, double y) {
    double FL = Structf_LLO(Q,xbj);
    double FT = Structf_TLO(Q,xbj);
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigmar = F2 - fy*FL;
    return sigmar;
}

// STRUCTURE FUNCTIONS -------------------------
double Structf_LLO ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LLOp( Q , xbj );
    return FL;
}

double Structf_TLO ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TLOp( Q , xbj );
    return FT;
}


///===========================================================================================
// CUBA WRAP
namespace cuba_config{
    int verbose=0;
    //maxeval=nlodis_config::CUBA_MAXEVAL;
    //double epsrel=nlodis_config::CUBA_EPSREL;
    double epsabs=0;
}

void Cuba(string method, int ndim, integrand_t integrand,
    void *userdata, double *integral, double *error, double *prob) {
    // common arguments
    int ncomp=1, nvec=1, seed=0, mineval=0, last=4;
    int nregions, neval, fail;
    void *spin=NULL;
    char *statefile=NULL;
    if(method=="vegas"){
    // Vegas-specific arguments
    int nstart=1000, nincrease=500, nbatch=1000, gridno=0;
    Vegas(ndim,ncomp,integrand,userdata,nvec,nlodis_config::CUBA_EPSREL,
        cuba_config::epsabs,cuba_config::verbose,seed,mineval,
        nlodis_config::CUBA_MAXEVAL,nstart,nincrease,nbatch,gridno,statefile,
        spin,&neval,&fail,integral,error,prob);
    }
    else if(method=="suave"){
    // Suave-specific arguments
    int nnew=1e3, nmin=2; // nnew=10e3
    double flatness=25; //25;
    Suave(ndim,ncomp,integrand,userdata,nvec,nlodis_config::CUBA_EPSREL,
        cuba_config::epsabs,cuba_config::verbose | last,seed,mineval,
        nlodis_config::CUBA_MAXEVAL,nnew,nmin,flatness,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
    else if(method=="divonne"){
    if(ndim==1) ndim=2;
    // Divonne-specific arguments
    int key1=1*47, key2=1, key3=1, maxpass=5, ngiven=0, nextra=0;
    double border=1e-8, maxchisq=10, mindeviation=0.25;
    Divonne(ndim,ncomp,integrand,userdata,nvec,nlodis_config::CUBA_EPSREL,
        cuba_config::epsabs,cuba_config::verbose,seed,mineval,
        nlodis_config::CUBA_MAXEVAL,key1,key2,key3,maxpass,border,maxchisq,
        mindeviation,ngiven,ndim,NULL,nextra,NULL,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
    else if(method=="cuhre"){
    if(ndim==1) ndim=2;
    // Cuhre-specific arguments
    int key=0;
    Cuhre(ndim,ncomp,integrand,userdata,nvec,nlodis_config::CUBA_EPSREL,
        cuba_config::epsabs,cuba_config::verbose | last,mineval,
        nlodis_config::CUBA_MAXEVAL,key,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
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


struct Userdata{
    double Q, xbj;
    double qMass;
    double icX0;
};

/*
// --- L L L --- LO -------- L L L --- LO -------- L L L --- LO -----------
*/
double ILLO(double Q, double z1, double x01sq) {
    double bessel_inner_fun = Q*sqrt(z1*(1.0-z1)*x01sq);
    double res = 0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        res = 0;
    }else{
        res = 4.0*Sq(Q)*Sq(z1)*Sq(1.0-z1)*Sq(gsl_sf_bessel_K0( bessel_inner_fun ));
    }   
    return res;
}

int integrand_ILLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=x01*x01;
    double Xrpdty_lo = Xrpdty_LO(xbj, Q*Q);

    double res=(1.0-(Sr(x01,Xrpdty_lo)))*(ILLO(Q,z1,x01sq))*x01;
        *f=res;
    return 0;
}

double LLOp(double Q, double x) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Cuba(cubaMethod,ndim,integrand_ILLOp,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}


/*
// --- T T T --- LO -------- T T T --- LO -------- T T T --- LO -----------
*/

double ITLO(double Q, double z1, double x01sq) {
    double bessel_inner_fun = Q*sqrt(z1*(1.0-z1)*x01sq);
    double res = 0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        res = 0;
    }else{
        res = Sq(Q)*z1*(1.0-z1)*(1.0-2.0*z1+2.0*Sq(z1))*Sq(gsl_sf_bessel_K1(bessel_inner_fun));
    }
    return res;
}

int integrand_ITLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=x01*x01;
    double Xrpdty_lo = Xrpdty_LO(xbj, Sq(Q));

    double res=(1.0-(Sr(x01,Xrpdty_lo)))*(ITLO(Q,z1,x01sq))*x01;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double TLOp(double Q, double x) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Cuba(cubaMethod,ndim,integrand_ITLOp,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}



/*
// Fitter tooling
*/
void InverseDipoleFitter::AddDataset(Data& d)
{
    datasets.push_back(&d);
}

InverseDipoleFitter::InverseDipoleFitter(MnUserParameters parameters_)
{
    parameters = parameters_;
    
    cubaMethod = "vegas";  // Default choise for Cuba
}


double InverseDipoleFitter::operator()(const std::vector<double>& par) const
{
    double chisqr = 0;
    // params to IC
    // double qs0sqr = 0.165;
    double alphas_scaling = 6.35;
    double anomalous_dimension  = 1.135;
    double icx0_nlo_impfac  = 1e-2;
    double icx0_bk  = 1e-2;
    double initialconditionY0  = 0;
    double icTypicalPartonVirtualityQ0sqr = 1;
    double old_sigma02 = 1;
    double qMass_charm = 1.27;
    double qMass_b = 4.750; // GeV, pole mass scheme val
    bool useMasses = nlodis_config::USE_MASSES;
    bool useCharm = false;
    
    std::vector< double > Svals = par;
    for(int i=0; i < Svals.size(); i++){
        auto& [r, x, S] = dipoleGrid[i];
        S = Svals[i];
    }

    cout << "=== Computing Reduced Cross sections ===" << endl;

    /*
     * Loop over datapoints and compute theoretical predictions
     */
    int points=0, totalpoints = 0;
    for (unsigned int dataset=0; dataset<datasets.size(); dataset++)
        totalpoints += datasets[dataset]->NumOfPoints();
    
    double fitsigma0 = 1;
    std::vector<double> datavals(totalpoints);
    std::vector<double> dataerrs(totalpoints);
    std::vector<double> thdata(totalpoints);  
    std::vector<double> var_xbj(totalpoints);
    std::vector<double> var_y(totalpoints);
    std::vector<double> var_qsqr(totalpoints);
    
    ///TODO: Does not yet fully support weight factors ,
    // after I have written a separate code that automatically finds optimal sigma02

    // These loops are trivially parallerizable
    // We only parallerize the inner loop where we have about
    // 250 points (total sigmar) and 50 points (charm)
    for (unsigned int dataset=0; dataset<datasets.size(); dataset++)
    {
#ifdef PARALLEL_CHISQR
        //reduction(+:chisqr)
    #pragma omp parallel for schedule(dynamic) reduction(+:points)
#endif
        for (int i=0; i<datasets[dataset]->NumOfPoints(); i++)
        {
            // Progress indication during fitting.
            #pragma omp critical
            cout << "\r" << i+1 << "/" << totalpoints << flush;

            // Index for this in the final data array
            int dataind=0;
            for (int dseti=0; dseti < dataset; dseti++)
                dataind += datasets[dseti]->NumOfPoints();
            dataind += i;
            
            
            
            double xbj      = datasets[dataset]->xbj(i);
            double y        = datasets[dataset]->y(i);              // inelasticity
            double Q2       = datasets[dataset]->Qsqr(i);
            double Q        = sqrt(Q2);
            double sigmar   = datasets[dataset]->ReducedCrossSection(i);
            double sigmar_err = datasets[dataset]->ReducedCrossSectionError(i);
            
            datavals[dataind] = sigmar;
            dataerrs[dataind] = sigmar_err;
            var_xbj[dataind] = xbj;
            var_y[dataind] = y;
            var_qsqr[dataind] = Q2;

            double theory=0, theory_charm=0;
            int calccount=0;
            if (!computeNLO && !useMasses) // Compute reduced cross section using leading order impact factors
            {
                theory = (fitsigma0)*SigmarLO(Q , xbj , y );
                ++calccount;
            }
            // if (!computeNLO && useMasses)
            // {
            //     theory=0;
            //     if (datasets[dataset]->OnlyCharm(i)==false)
            //     // if (true)
            //     {
            //         theory = (fitsigma0)*SigmaComputer.SigmarLOmass(Q , xbj , y, false );
            //     }
            //     if (xbj*(1.0 + 4.0*1.35*1.35/(Q*Q)) < 0.01 and useCharm)
            //     {
            //         theory_charm = (fitsigma0)*SigmaComputer.SigmarLOmass(Q , xbj*(1.0 + 4.0*1.35*1.35/(Q*Q)) , y, true ); // charm
            //         if (datasets[dataset]->OnlyCharm(i) == true)
            //             theory = theory_charm;
            //         else
            //             theory = theory + theory_charm;
            //     }
            //       ++calccount;
            // }
            
            // if (computeNLO && !UseSub) // UNSUB SCHEME Full NLO impact factors for reduced cross section
            // {
            //     if (useMasses){
            //         if (useBoundLoop){
            //             cout << "Z_2 bound dipole term not implemented yet with quark masses. EXIT." << endl;
            //             exit(1);
            //             // theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub_UniformZ2Bound(Q , xbj , y );
            //             ++calccount;}
            //         if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
            //             // cout << "This should be in use???" << "charm mass is: " << qMass_charm << endl;
            //             if (nlodis_config::MASS_SCHEME == nlodis_config::CHARM_ONLY){
            //                 theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub_massive(Q , xbj , y, qMass_charm );
            //             } else if (nlodis_config::MASS_SCHEME == nlodis_config::BEAUTY_ONLY){
            //                 theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub_massive(Q , xbj , y, qMass_b_var );
            //             } else if (nlodis_config::MASS_SCHEME == nlodis_config::LIGHT_PLUS_CHARM){
            //                 theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub(Q , xbj , y )
            //                          + (fitsigma0)*SigmaComputer.SigmarNLOunsub_massive(Q , xbj , y, qMass_charm );
            //             } else if (nlodis_config::MASS_SCHEME == nlodis_config::LIGHT_PLUS_CHARM_AND_BEAUTY){
            //                 theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub(Q , xbj , y )
            //                          + (fitsigma0)*SigmaComputer.SigmarNLOunsub_massive(Q , xbj , y, qMass_charm )
            //                          + (fitsigma0)*SigmaComputer.SigmarNLOunsub_massive(Q , xbj , y, qMass_b_var );
            //             }
            //             ++calccount;}
            //     }
            //     if (!useMasses){
            //         if (useBoundLoop){
            //             theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub_UniformZ2Bound(Q , xbj , y );
            //             ++calccount;}
            //         if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
            //             theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub(Q , xbj , y );
            //             //theory = (fitsigma0)*SigmaComputer.SigmarNLOsubRisto(Q , xbj , y );
            //             ++calccount;}
            //         if (UseSigma3){
            //             theory += (fitsigma0)*SigmaComputer.SigmarNLOunsub_sigma3(Q , xbj , y );
            //             }
            //     }
            // }
            
            // if (computeNLO && UseSub) // SUB SCHEME Full NLO impact factors for reduced cross section
            // {
            //     if (useMasses){
            //         cout << "SUB SCHEME NOT IMPLEMENTED WITH QUARK MASSES. EXIT." << endl;
            //         exit(1);
            //     }
            //     if (useBoundLoop){
            //         theory = (fitsigma0)*SigmaComputer.SigmarNLOsub_UniformZ2Bound(Q , xbj , y );
            //         ++calccount;}
            //     if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
            //         theory = (fitsigma0)*SigmaComputer.SigmarNLOsub(Q , xbj , y );
            //         //theory = (fitsigma0)*SigmaComputer.SigmarNLOsubRisto(Q , xbj , y );
            //         ++calccount;}
            // }
            
            thdata[dataind] = theory;


            if (calccount>1)
            {
              cerr << "ERROR: Multiple computations. abort." << "count="<< calccount << endl;
              theory = 99999999;
              thdata[dataind] = 99999999;
              exit(1);
            }

            if (std::isnan(theory) or std::isinf(theory))
            {
                cerr << "Warning: theory result " << theory << " with parameters " << PrintVector(par) << endl;
                theory = 99999999;
            }

            //chisqr += datasets[dataset]->Weight()*SQR( (theory+theory_charm - sigmar) / sigmar_err );
            points = points + datasets[dataset]->Weight();
        }
    }
    cout << endl;
    // Minimize sigma02
    std::vector<double> sigma02fit = FindOptimalSigma02(datavals,dataerrs, thdata);
    double chisqr_over_n = sigma02fit[1];
    double sigma02 = sigma02fit[0];
    // comparing to old sigma02
    // double old_chisq_over_n = minimiser_helper_chisqr_vec(old_sigma02, datavals, dataerrs, thdata);
    // Output for plotting
    if(false){
        for(int i=0; i<var_xbj.size(); i++){
        #pragma omp critical
        cout    << setw(10) << var_xbj[i]      << " "
                << setw(10) << var_qsqr[i]     << " "
                << setw(10) << var_y[i]        << " "
                << setw(10) << datavals[i]     << " "
                << setw(10) << dataerrs[i]     << " "
                << setw(10) << sigma02*thdata[i] /* << " "
                << setw(10) << sigma02*theory_charm <<*/ 
                << endl;
            }
    }
    cout    << endl 
            << "# Calculated chi^2/N = " << chisqr_over_n
            << " (N=" << points
            << "), parameters (" << PrintVector(par)
            << ", sigma02=" << sigma02
            << ")" << endl<<endl;
    // cout    << endl
    //         << "# Calculated chi^2/N = " << old_chisq_over_n
    //         << " (N=" << points
    //         << "), parameters (" << PrintVector(par)
    //         << ", sigma02=" << old_sigma02
    //         << ")" << endl<<endl;
    

    // delete DipolePointer;
    return chisqr_over_n * points;  // Return chisqr for Minuit error estimation.
};





int main( int argc, char* argv[] )
{
    gsl_set_error_handler(&ErrHandlerCustom2);

    // NLO DIS SIGMA_R COMPUTATION CONFIGS
    nlodis_config::USE_MASSES = true;

    nlodis_config::CUBA_EPSREL = 1e-3;
    nlodis_config::CUBA_MAXEVAL= 1e7;
    nlodis_config::MINR = 1e-6;
    nlodis_config::MAXR = 30;
    nlodis_config::PRINTDATA = true;
    bool useNLO = false;
    bool computeNLO = useNLO;
    string cubaMethod = "vegas";
    // string cubaMethod = "suave";

    config::NO_K2 = true;  // Do not include numerically demanding full NLO part
    config::KINEMATICAL_CONSTRAINT = config::KC_NONE;

    config::VERBOSE = true;
    //config::RINTPOINTS = 512/4;
    //config::THETAINTPOINTS = 512/4;

    config::INTACCURACY = 5e-3;//0.02;
    // config::INTACCURACY = 5e-3; // highacc def1
    // config::INTACCURACY = 20e-3; // quick low acc
    config::MINR = 1e-6;
    config::MAXR = 30;
    // config::MINR = 1e-4;  // faster lower accuracy limits
    // config::MAXR = 20;
    config::RPOINTS = 100;
    config::DE_SOLVER_STEP = 0.4; // Rungekutta step

    // Constants
    config::NF=3;   // Only light quarks
    config::LAMBDAQCD = 0.241;


    Data data;
    data.SetMinQsqr(0.75);
    data.SetMaxQsqr(50);
    data.SetMaxX(0.01);


    MnUserParameters parameters;

    bool useSUB, useResumBK, useKCBK, useImprovedZ2Bound, useBoundLoop;
    bool useSigma3 = false;
    string helpstring = "Argument order: SCHEME BK RC useImprovedZ2Bound useBoundLoop [Qs0 C^2 gamma] X0_if X0_bk e_c Q0sq Y0 eta0\nsub/unsub/unsub+ resumbk/trbk/lobk parentrc/guillaumerc/fixedrc z2improved/z2simple z2boundloop/unboundloop";
    string string_sub, string_bk, string_rc;
    if (argc<2){ cout << helpstring << endl; return 0;}
    // Argv[0] is the name of the program

    nlodis_config::MASS_SCHEME = nlodis_config::MASSLESS;
    string_sub = string(argv [1]);
    if (string(argv [1]) == "sub"){
        useSUB = true;
        nlodis_config::USE_MASSES = false;
        nlodis_config::SUB_SCHEME = nlodis_config::SUBTRACTED;
    } else if (string(argv [1]) == "unsub"){
        nlodis_config::SUB_SCHEME = nlodis_config::UNSUBTRACTED;
        useSUB = false;
        useSigma3 = false;
        nlodis_config::USE_MASSES = false;
    } else if (string(argv [1]) == "uncc"){
        nlodis_config::SUB_SCHEME = nlodis_config::UNSUBTRACTED;
        useSUB = false;
        useSigma3 = false;
        nlodis_config::MASS_SCHEME = nlodis_config::CHARM_ONLY;
    } else if (string(argv [1]) == "uncc2"){
        nlodis_config::SUB_SCHEME = nlodis_config::UNSUBTRACTED;
        useSUB = false;
        useSigma3 = false;
        nlodis_config::MASS_SCHEME = nlodis_config::CHARM_ONLY;
        nlodis_config::PERF_MODE = nlodis_config::MASSIVE_EXPLICIT_BESSEL_DIM_REDUCTION;
    } else if (string(argv [1]) == "unbb"){
        nlodis_config::SUB_SCHEME = nlodis_config::UNSUBTRACTED;
        useSUB = false;
        useSigma3 = false;
        nlodis_config::MASS_SCHEME = nlodis_config::BEAUTY_ONLY;
    } else if (string(argv [1]) == "unbb2"){
        nlodis_config::SUB_SCHEME = nlodis_config::UNSUBTRACTED;
        useSUB = false;
        useSigma3 = false;
        nlodis_config::MASS_SCHEME = nlodis_config::BEAUTY_ONLY;
        nlodis_config::PERF_MODE = nlodis_config::MASSIVE_EXPLICIT_BESSEL_DIM_REDUCTION;
    } else if (string(argv [1]) == "unlpc"){
        nlodis_config::SUB_SCHEME = nlodis_config::UNSUBTRACTED;
        useSUB = false;
        useSigma3 = false;
        nlodis_config::MASS_SCHEME = nlodis_config::LIGHT_PLUS_CHARM;
    } else if (string(argv [1]) == "unlpcb"){
        nlodis_config::SUB_SCHEME = nlodis_config::UNSUBTRACTED;
        useSUB = false;
        useSigma3 = false;
        nlodis_config::MASS_SCHEME = nlodis_config::LIGHT_PLUS_CHARM_AND_BEAUTY;
    } else if (string(argv [1]) == "unlpcb2"){
        nlodis_config::SUB_SCHEME = nlodis_config::UNSUBTRACTED;
        useSUB = false;
        useSigma3 = false;
        nlodis_config::MASS_SCHEME = nlodis_config::LIGHT_PLUS_CHARM_AND_BEAUTY;
        nlodis_config::PERF_MODE = nlodis_config::MASSIVE_EXPLICIT_BESSEL_DIM_REDUCTION;
    } else {cout << helpstring << endl; return -1;}

    if (nlodis_config::MASS_SCHEME == nlodis_config::CHARM_ONLY){
        // data.LoadData("./data/hera_combined_sigmar_cc.txt", TOTAL); // old charm data
        data.LoadData("./data/hera_II_combined_sigmar_cc.txt", TOTAL); // newer charm data
    } else if (nlodis_config::MASS_SCHEME == nlodis_config::BEAUTY_ONLY){
        data.LoadData("./data/hera_II_combined_sigmar_b.txt", TOTAL); // newer bottom data
    } else {
        data.LoadData("./data/hera_combined_sigmar.txt", TOTAL); // older total data
        // data.LoadData("./data/hera_II_combined_sigmar.txt", TOTAL); // newer total data
    }

    string_bk = string(argv [2]);
    if (string(argv [2]) == "resumbk"){
            config::EULER_METHOD = false;    // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = true;       // Resum doulbe logs
            config::RESUM_SINGLE_LOG = true; // Resum single logs
            config::KSUB = 0.65;             // Optimal value for K_sub
            config::NO_K2 = true;            // Do not include numerically demanding full NLO part
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_RESUM;
    }else if (string(argv [2]) == "kcbk"){
            config::EULER_METHOD = true;     // Kinematical constraint requires this
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_BEUF_K_PLUS;
            config::DE_SOLVER_STEP = 0.05;  //0.02; // Euler method requires smaller step than RungeKutta!
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_KCBK_BEUF;
    }else if (string(argv [2]) == "trbk"){  // Target Rapidity BK
            config::EULER_METHOD = true;     // Kinematical constraint requires this
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_EDMOND_K_MINUS;
            config::DE_SOLVER_STEP = 0.05;  //0.02; // Euler method requires smaller step than RungeKutta!
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_TRBK_EDMOND;
            nlodis_config::TRBK_RHO_PRESC = nlodis_config::TRBK_RHO_RQ0;
            // nlodis_config::TRBK_RHO_PRESC = nlodis_config::TRBK_RHO_QQ0;
    }else if (string(argv [2]) == "nlobk"){
            config::EULER_METHOD = false;    // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = true;       // Resum doulbe logs
            config::RESUM_SINGLE_LOG = true; // Resum single logs
            config::KSUB = 0.65;             // Optimal value for K_sub
            config::NO_K2 = false;            // Do not include numerically demanding full NLO part
    }else if (string(argv [2]) == "lobk"){
            config::EULER_METHOD = false;   // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_NONE;
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_LOBK_EXPLICIT;
    }else if (string(argv [2]) == "lobkold"){
            config::EULER_METHOD = false;   // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_NONE;
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_LOBK_Z2TOZERO;
    } else {cout << helpstring << endl; return -1;}

    string_rc = string(argv [3]);
    if (string(argv [3]) == "parentrc" or string(argv [3]) == "pdrc"){
            config::RC_LO = config::PARENT_LO;
            config::RC_NLO = config::PARENT_NLO;
            config::RESUM_RC = config::RESUM_RC_PARENT;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_PARENT;
    } else if (string(argv [3]) == "guillaumerc" or string(argv [3]) == "gbrc"){
            config::RC_LO = config::GUILLAUME_LO;
            config::RESUM_RC = config::RESUM_RC_GUILLAUME;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_GUILLAUME;
    } else if (string(argv [3]) == "fixedrc" or string(argv [3]) == "fc"){
            config::RC_LO = config::FIXED_LO;
            config::RESUM_RC = config::RESUM_RC_FIXED;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_FIXED;
    } else if (string(argv [3]) == "smallestrc" or string(argv [3]) == "sdrc"){
            config::RC_LO = config::SMALLEST_LO;
            config::RC_NLO = config::SMALLEST_NLO;
            config::RESUM_RC = config::RESUM_RC_SMALLEST;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_SMALLEST;
    } else if (string(argv [3]) == "balitskysmallrc" or string(argv [3]) == "balsdrc"){
            // Even though this coupling should be more realistic than smallest dipole alone,
            // the subtraction between the LO and qg terms is not exact. This shortcoming makes
            // this coupling less than ideal.
            config::RC_LO = config::BALITSKY_LO;
            config::RC_NLO = config::SMALLEST_NLO;
            config::RESUM_RC = config::RESUM_RC_SMALLEST;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_SMALLEST;
    }else {cout << helpstring << endl; return -1;}

    if (string(argv [4]) == "z2improved" or string(argv [4]) == "z2imp"){
        nlodis_config::Z2MINIMUM = nlodis_config::Z2IMPROVED;
        useImprovedZ2Bound = true;
    } else if (string(argv [4]) == "z2simple" or string(argv [4]) == "z2sim"){
        nlodis_config::Z2MINIMUM = nlodis_config::Z2SIMPLE;
        useImprovedZ2Bound = false;
    } else {cout << helpstring << endl; return -1;}

    if (string(argv [5]) == "z2boundloop" or string(argv [5]) == "z2b"){
        useBoundLoop = true;
    } else if (string(argv [5]) == "unboundloop" or string(argv [5]) == "unb"){
        useBoundLoop = false;
    } else {cout << helpstring << endl; return -1;}

    // loading dipole function initial shape
    AmplitudeLib* DipoleAmplitude_ptr; // Forward declaration of the dipole object to be initialized from a file or solved data.
    // string dipole_basename = "./out/dipoles/dipole";
    // string dipole_filename = dipole_basename
    //                          + "_" + string_bk
    //                          + "_" + string_rc
    //                          + "_x0bk" + std::to_string(icx0_bk)
    //                          + "_qs0sqr" + std::to_string(qs0sqr)
    //                          + "_asC^2" + std::to_string(alphas_scaling)
    //                          + "_gamma" + std::to_string(anomalous_dimension)
    //                          + "_ec" + std::to_string(e_c)
    //                          + "_eta0" + std::to_string(eta0)
    //                          + "_maxy" + std::to_string(maxy)
    //                          + "_euler" + std::to_string(config::EULER_METHOD)
    //                          + "_step" + std::to_string(config::DE_SOLVER_STEP)
    //                          + "_rpoints" + std::to_string(config::RPOINTS)
    //                          + "_rminmax" + std::to_string(config::MINR) + "--" + std::to_string(config::MAXR)
    //                          + "_intacc" + std::to_string(config::INTACCURACY) ;
    string dipole_filename = string(argv [6]);
    
    // generate publishable filenames
    // string bk_name = (string_bk == "trbk") ? "tbk" : string_bk;
    // string rc_name = (string_rc == "sdrc") ? "bal+sd" : "parent";
    // string Y0_valu = (icx0_bk == 1.0) ? "0.00" : std::to_string((int)(std::log(1./icx0_bk) * 100 + .5) / 100.0);
    // Y0_valu.erase(Y0_valu.find_last_of(".") + 3, std::string::npos);
    // string dipole_filename = dipole_basename
    //                          + "-" + bk_name
    //                          + "-" + dataname
    //                          + "-" + rc_name
    //                          + "-" + Y0_valu
    //                          + ".dip";
    if (FILE *file = fopen(dipole_filename.c_str(), "r")) {
        cout << "# Previously saved dipole file found: " << dipole_filename << endl;
        DipoleAmplitude_ptr = new AmplitudeLib(dipole_filename);      // read data from existing file.
        fclose(file);
    } 
    // else {
    //     solver.Solve(maxy);     // Solve up to maxy since specified dipole datafile was not found.
    //     solver.GetDipole()->Save(dipole_filename);
    //     cout << "# Saved dipole to file: "<< dipole_filename << endl;
    //     DipoleAmplitude_ptr = new AmplitudeLib(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());
    // }   
    AmplitudeLib DipoleAmplitude(*DipoleAmplitude_ptr);
    DipoleAmplitude.SetInterpolationMethod(LINEAR_LINEAR);
    double icx0_bk = 1e-2;
    DipoleAmplitude.SetX0(icx0_bk);
    DipoleAmplitude.SetOutOfRangeErrors(false);
    // AmplitudeLib *DipolePointer = &DipoleAmplitude;
    invdip_config::DipolePointer = &DipoleAmplitude; // replace the above with a global pointer
    // Discretizing dipole amplitude
    int ngrid = 2; // r * Y grid size
    double rmin = config::MINR;
    double rmax = config::MAXR;
    double xmin = 1e-5; //1.85E-05
    double xmax = icx0_bk; //bk ic
    
    // INITIALIZE DISCRETE DIPOLE
    // double rstep = (rmax/rmin)/ngrid;
    // double rstep = std::pow(rmax/rmin,1./((double)ngrid));
    double rstep = (rmax-rmin)/((double)ngrid);
    double xstep = std::pow(xmax/xmin,1./((double)ngrid));
    // generate grid
    double r, x;
    r = rmin;
    x = xmax;
    for (int i=0; i<ngrid+1; i++){
        for (int j=0; j<ngrid+1; j++){
            // r = rmin*std::pow(rstep,(double)i);
            r = rmin+rstep*(double)i;
            x = xmax/std::pow(xstep,(double)j);
            dipoleGrid.emplace_back(r,x,DipolePointer->S(r,x));
            rvals.push_back(r);
            xvals.push_back(x);
            Svals.push_back(DipolePointer->S(r,x));
            cout << "r " << r << " x " << x << " Srx " << DipolePointer->S(r,x) << endl; 
        }
    }
    // exit(0);

    double old_sigma02 = 1.;
    double mass_charm = 1.27; // MSbar value as default ----- // OLD WAS 1.35;
    double mass_bottom = 4.75; // default pole scheme value used before fitting.


    for (auto elem : dipoleGrid)
    {
        auto [r, x, S] = elem;
        string param = "S_r" + std::to_string(r) + "x_" + std::to_string(x);
        if (r == rmin){
	        parameters.Add(param, S);
        } else {
	    parameters.Add(param, S, 0.05);
        parameters.SetLimits(param, 0.0, 1.0);
        }
    }
    

    InverseDipoleFitter fitter(parameters);
    fitter.AddDataset(data);
    fitter.AddDipoleGrid(dipoleGrid);
    fitter.SetNLO(useNLO);
    fitter.SetSUB(useSUB);
    fitter.SetSigma3(useSigma3);
    fitter.UseImprovedZ2Bound(useImprovedZ2Bound);
    fitter.UseConsistentlyBoundLoopTerm(useBoundLoop);
    fitter.SetCubaMethod(cubaMethod);

    cout << std::boolalpha;
    cout    << "# === Perturbative settings ===" << endl
            << "# Use masses: " << nlodis_config::USE_MASSES << ", scheme:" << nlodis_config::MASS_SCHEME << endl
            << "# Settings: " << string_sub << " (scheme), " << string_bk << ", " << string_rc << endl
            << "# Use LOBK (DL,SL==false): " << (!(config::RESUM_DLOG) 
                                    and !(config::RESUM_SINGLE_LOG)) << endl
            << "# Use ResumBK (DL,SL==true,KC_NONE): " << ((config::RESUM_DLOG) 
                                    and (config::RESUM_SINGLE_LOG)
                                    and (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)
                                    and (config::NO_K2 == true)) << endl
            << "# Use NLOBK (DL,SL==true,KC_NONE,NO_K2==false): " << ((config::RESUM_DLOG) 
                                    and (config::RESUM_SINGLE_LOG)
                                    and (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)
                                    and (config::NO_K2 == false)) << endl
            << "# KinematicalConstraint / target eta0 BK: " << config::KINEMATICAL_CONSTRAINT << " (0 BEUF_K_PLUS, 1 EDMOND_K_MINUS, 2 NONE)" << endl
            << "# Target eta0 RHO shift: " << nlodis_config::TRBK_RHO_PRESC << " (0 TRBK_RHO_DISABLED, 1 TRBK_RHO_QQ0, 2 TRBK_RHO_RQ0)" << endl
            << "# Running Coupling: (RC_LO):    " << config::RC_LO << " (0 fc, 1 parent, 2 parent_beta, 3 smallest, 4 balitsky, 5 frac, 6 guillaume)" << endl
            << "# Running Coupling: (RC_NLO):    " << config::RC_NLO << " (0 fc, 1 parent, 2 smallest)" << endl
            << "# Running Coupling: (RESUM_RC): " << config::RESUM_RC << " (0 fc, 1 balitsky, 2 parent, 3 smallest, 4 guillaume)" << endl
            << "# Running Coupling: (RC_DIS):   " << nlodis_config::RC_DIS << " (0 fc, 1 parent, 2 smallest, 3 guillaume)" << endl
            << "# Use NLOimpact: " << useNLO << endl
            << "# Use SUBscheme: " << useSUB << endl
            << "# Use Sigma3: " << useSigma3 << endl
            << "# Use improved Z2 bound: " << useImprovedZ2Bound << endl
            << "# Use Z2 loop term: " << useBoundLoop << endl
            << "# Cuba MC: " << cubaMethod
                << ", Cuba eps = " << nlodis_config::CUBA_EPSREL
                << ", Cuba maxeval = " << (float)nlodis_config::CUBA_MAXEVAL
                << ", Cuba perf scheme = " << nlodis_config::PERF_MODE
                << endl
            << "# config::INTACCURACY = " << config::INTACCURACY
                << ", config::RPOINTS = " << config::RPOINTS
                << ", config::DE_SOLVER_STEP = " << config::DE_SOLVER_STEP
                << ", config::{MINR, MAXR}, nlodis_config::{MINR, MAXR} = " << config::MINR << " " << config::MAXR << " " << nlodis_config::MINR << " " << nlodis_config::MAXR
                << endl;
    cout << "=== Initial parameters ===" << endl;
    cout << parameters << endl;
    if(nlodis_config::VERBOSE) cout << "=== Starting fit ===" << endl;

    MnMigrad fit(fitter, parameters, 0);
    // MnMinimize fit(fitter, parameters);
    fit.SetPrecision(1e-6);
    FunctionMinimum min = fit();
    std::cout<<"minimum: "<<min<<std::endl;

    cout << std::boolalpha;
    cout    << "# === Perturbative settings ===" << endl
            << "# Use masses: " << nlodis_config::USE_MASSES << ", scheme:" << nlodis_config::MASS_SCHEME << endl
            << "# Settings: " << string_sub << " (scheme), " << string_bk << ", " << string_rc << endl
            << "# Use LOBK (DL,SL==false): " << (!(config::RESUM_DLOG) 
                                    and !(config::RESUM_SINGLE_LOG)) << endl
            << "# Use ResumBK (DL,SL==true,KC_NONE): " << ((config::RESUM_DLOG) 
                                    and (config::RESUM_SINGLE_LOG)
                                    and (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)
                                    and (config::NO_K2 == true)) << endl
            << "# Use NLOBK (DL,SL==true,KC_NONE,NO_K2==false): " << ((config::RESUM_DLOG) 
                                    and (config::RESUM_SINGLE_LOG)
                                    and (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)
                                    and (config::NO_K2 == false)) << endl
            << "# KinematicalConstraint / target eta0 BK: " << config::KINEMATICAL_CONSTRAINT << " (0 BEUF_K_PLUS, 1 EDMOND_K_MINUS, 2 NONE)" << endl
            << "# Target eta0 RHO shift: " << nlodis_config::TRBK_RHO_PRESC << " (0 TRBK_RHO_DISABLED, 1 TRBK_RHO_QQ0, 2 TRBK_RHO_RQ0)" << endl
            << "# Running Coupling: (RC_LO):    " << config::RC_LO << " (0 fc, 1 parent, 2 parent_beta, 3 smallest, 4 balitsky, 5 frac, 6 guillaume)" << endl
            << "# Running Coupling: (RC_NLO):    " << config::RC_NLO << " (0 fc, 1 parent, 2 smallest)" << endl
            << "# Running Coupling: (RESUM_RC): " << config::RESUM_RC << " (0 fc, 1 balitsky, 2 parent, 3 smallest, 4 guillaume)" << endl
            << "# Running Coupling: (RC_DIS):   " << nlodis_config::RC_DIS << " (0 fc, 1 parent, 2 smallest, 3 guillaume)" << endl
            << "# Use NLOimpact: " << useNLO << endl
            << "# Use SUBscheme: " << useSUB << endl
            << "# Use Sigma3: " << useSigma3 << endl
            << "# Use improved Z2 bound: " << useImprovedZ2Bound << endl
            << "# Use Z2 loop term: " << useBoundLoop << endl
            << "# Cuba MC: " << cubaMethod
                << ", Cuba eps = " << nlodis_config::CUBA_EPSREL
                << ", Cuba maxeval = " << (float)nlodis_config::CUBA_MAXEVAL
                << ", Cuba perf scheme = " << nlodis_config::PERF_MODE
                << endl
            << "# config::INTACCURACY = " << config::INTACCURACY
                << ", config::RPOINTS = " << config::RPOINTS
                << ", config::DE_SOLVER_STEP = " << config::DE_SOLVER_STEP
                << ", config::{MINR, MAXR}, nlodis_config::{MINR, MAXR} = " << config::MINR << " " << config::MAXR << " " << nlodis_config::MINR << " " << nlodis_config::MAXR
                << endl;


    return 0;
}
