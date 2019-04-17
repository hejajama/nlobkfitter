#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <csignal>
#include <ctime>
#include <gsl/gsl_errno.h>
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
#include "nlodis_config.hpp"
#include "nlodissigmar.hpp"

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnScan.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>

 #include <gsl/gsl_errno.h>

 #include "gitsha1.h"

using namespace std;
using namespace ROOT::Minuit2;
//void ErrHandler(const char * reason,const char * file,int line,int gsl_errno);
int errors_count;
void ErrHandlerCustom(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{   // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    // 15: underflows

    if (gsl_errno == 11) return; // ignore max subdivision errors
    if (gsl_errno == 18) return; // roundoff errors from small r?
    // if (gsl_errno == 15 or gsl_errno == 16) return; // underfows and overflows?
    // Ugly hack, comes from the edges of the z integral in virtual_photon.cpp
    // Overflows come from IPsat::bint when it is done analytically
    // Hope is that these errors are handled correctly everywhere
    errors_count++;
    std::cerr << file << ":"<< line <<": Error " << errors_count << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}

int main( int argc, char* argv[] )
{
    gsl_set_error_handler(&ErrHandlerCustom);

        // NLO DIS SIGMA_R COMPUTATION CONFIGS
        nlodis_config::CUBA_EPSREL = 10e-3;
        nlodis_config::CUBA_MAXEVAL= 1e7;
        nlodis_config::MINR = 1e-5;
        nlodis_config::MAXR = 50;
        nlodis_config::PRINTDATA = true;
        bool useNLO = true;
        //bool useSUB = true;               // Set by a command line switch in swarmscan
        //bool useImprovedZ2Bound = true;   // Set by a command line switch in swarmscan
        //bool useBoundLoop = true;         // Set by a command line switch in swarmscan
        string cubaMethod = "suave";
        //bool useResumBK = true;
        //bool useKCBK = false;
        //if (useResumBK == true and useKCBK == true) {cout << "Both ResumBK and KCBK enabled, exitting." << endl; return -1;}

        // Constants
        config::NF=3;   // Only light quarks
        config::LAMBDAQCD = 0.241;
        
        config::VERBOSE = true;
        config::RINTPOINTS = 512;
        config::THETAINTPOINTS = 512;

        nlodis_config::RC_DIS = nlodis_config::DIS_RC_GUILLAUME;
        config::RC_LO = config::GUILLAUME_LO;// FIXED_,PARENT_,PARENT_BETA_,SMALLEST_,BALITSKY_,FRAC_,GUILLAUME_,
        config::RESUM_RC = config::RESUM_RC_GUILLAUME; // _BALITSKY,_PARENT,_SMALLEST,_GUILLAUME,
        config::RESUM_DLOG = true; // Resum doulbe logs
        config::RESUM_SINGLE_LOG = true; // Resum single logs
        config::LO_BK = true;  // Solve LO BK with running coupling, overrides RESUM settings
        config::KSUB = 0.65;  // Optimal value for K_sub
        config::NO_K2 = true;  // Do not include numerically demanding full NLO part
        config::INTACCURACY = 10e-3;//0.02;
        config::MCINTACCURACY = 10e-3;//0.02;
        config::MCINTPOINTS = 1e7;
        config::MINR = 1e-5;
        config::MAXR = 50;
        config::RPOINTS = 100;
        config::DE_SOLVER_STEP = 0.4; // Rungekutta step
        // config::DE_SOLVER_STEP = 0.05; // Euler step


    // READING RUN CONFIGURATION FROM THE STDIN
    bool useSUB, useResumBK, useKCBK, useImprovedZ2Bound, useBoundLoop, useSigma3;
    string helpstring = "Argument order: SCHEME BK RC useImprovedZ2Bound useBoundLoop Q C^2 X0 gamma Q0sq Y0 eta0\nsub/unsub/unsub+ resumbk/trbk/lobk parentrc/guillaumerc/fixedrc z2improved/z2simple z2boundloop/unboundloop";
    string string_sub, string_bk, string_rc;
    if (argc<2){ cout << helpstring << endl; return 0;}
    // Argv[0] is the name of the program

    string_sub = string(argv [1]);
    if (string(argv [1]) == "sub"){
      useSUB = true;
    } else if (string(argv [1]) == "unsub"){
      useSUB = false;
      useSigma3 = false;
    } else if (string(argv [1]) == "unsub+"){
      useSUB = false;
      useSigma3 = true;
    } else {cout << helpstring << endl; return -1;}

    string_bk = string(argv [2]);
    if (string(argv [2]) == "resumbk"){
            config::LO_BK = false;  // Solve LO BK with running coupling, overrides RESUM settings
    }else if (string(argv [2]) == "kcbk"){
            cout << endl << "KCBK (old) unsupported at this time!" << endl << endl;
            exit(1);
            config::LO_BK = true;  // Solve (kinematic / delay) LO BK with running coupling, overrides RESUM settings
            config::EULER_METHOD            = true;        // Kinematical constraint requires this
            config::KINEMATICAL_CONSTRAINT  = true;
            config::DE_SOLVER_STEP = 0.08; //0.02; // Euler method requires smaller step than RungeKutta!
    }else if (string(argv [2]) == "trbk"){  // Target Rapidity BK
            config::LO_BK = true;  // Solve (kinematic / delay) LO BK with running coupling, overrides RESUM settings
            config::EULER_METHOD            = true;        // Kinematical constraint requires this
            config::KINEMATICAL_CONSTRAINT  = true;
            config::DE_SOLVER_STEP = 0.08; //0.02; // Euler method requires smaller step than RungeKutta!
    }else if (string(argv [2]) == "lobk"){
            config::LO_BK = true;
    } else {cout << helpstring << endl; return -1;}

    string_rc = string(argv [3]);
    if (string(argv [3]) == "parentrc"){
            config::RC_LO = config::PARENT_LO;
            config::RESUM_RC = config::RESUM_RC_PARENT;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_PARENT;
    } else if (string(argv [3]) == "guillaumerc"){
            config::RC_LO = config::GUILLAUME_LO;
            config::RESUM_RC = config::RESUM_RC_GUILLAUME;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_GUILLAUME;
    } else if (string(argv [3]) == "fixedrc"){
            if (config::LO_BK == false){
                cout << "There is no fixed coupling option for RESUM BK" << endl;
                exit(1); }
            config::RC_LO = config::FIXED_LO;
            // config::RESUM_RC = config::RESUM_RC_FIXED; // no such option
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_FIXED;
    } else {cout << helpstring << endl; return -1;}

    if (string(argv [4]) == "z2improved"){
      useImprovedZ2Bound = true;
    } else if (string(argv [4]) == "z2simple"){
      useImprovedZ2Bound = false;
    } else {cout << helpstring << endl; return -1;}

    if (string(argv [5]) == "z2boundloop"){
      useBoundLoop = true;
    } else if (string(argv [5]) == "unboundloop"){
      useBoundLoop = false;
    } else {cout << helpstring << endl; return -1;}

    int argi=5; argi++;
    // reading fit initial condition parameters:
    double icqs0sq   = stod( argv [argi] ); argi++;
    // double icsigma0  = stod( argv [argi] ); argi++;
    double iccsq     = stod( argv [argi] ); argi++;
    double icx0      = stod( argv [argi] ); argi++;
    double icgamma   = stod( argv [argi] ); argi++;
    double icQ0sq    = stod( argv [argi] ); argi++;
    double icY0      = stod( argv [argi] ); argi++;
    double icEta0    = stod( argv [argi] ); argi++;


    MnUserParameters parameters;
    // Fit parameters, first value is starting value, second is uncertainty
        //if(argc = 13){
        parameters.Add("qs0sqr",                icqs0sq , 0.2);
        // parameters.Add("qs0sqr",                icqs0sq);
        // parameters.Add("fitsigma0",             icsigma0, 10.0); // 1mb = 2.568 GeV² // (2.568)*16.36
        parameters.Add("alphascalingC2",        iccsq,    20.0);
        // parameters.Add("alphascalingC2",        iccsq);
        parameters.Add("e_c", 1.0 );
        parameters.Add("anomalous_dimension",   icgamma);
        parameters.Add("initialconditionX0",    icx0 );
        parameters.Add("initialconditionY0",    icY0 );
        parameters.Add("icTypicalPartonVirtualityQ0sqr", icQ0sq );
        parameters.Add("eta0", icEta0 );
        //}
        /*else{
        cout << "INSUFFICIENT RUN PARAMETERS, running a preset test run:" << endl;
        useSUB = true;
        useSigma3 = false;            
        config::LO_BK = false;  // Solve RESUM BK with running coupling
        useImprovedZ2Bound = true;  
        useBoundLoop = true;        
        parameters.Add("qs0sqr", 0.1 );
        parameters.Add("fitsigma0", 20.0 ); // 1mb = 2.568 GeV² // (2.568)*16.36
        parameters.Add("alphascalingC2", 0.315 );
        parameters.Add("e_c", 1.0 );
        parameters.Add("anomalous_dimension", 1.0 );
        parameters.Add("initialconditionX0", 1.0 );
        parameters.Add("initialconditionY0", 0. );
        parameters.Add("icTypicalPartonVirtualityQ0sqr", 1.0 ); }
        */

        //
        // Set limits
        //
        parameters.SetLimits("qs0sqr",              0.001 , 0.7);
        // parameters.SetLimits("fitsigma0",           1.0   , 40.0);
        parameters.SetLimits("alphascalingC2",      0.1  ,	30.0);
        //parameters.SetLimits("e_c",                 0.4   , 10.0);
        //parameters.SetLimits("anomalous_dimension", 0.1   ,	2.0);
        //parameters.SetLimits("initialconditionX0",  0.01  ,	10.0);
        //parameters.SetLimits("initialconditionY0",  0.01  ,	10.0);


    Data data;
    data.SetMinQsqr(0.75);
    data.SetMaxQsqr(50);
    data.SetMaxX(0.01);
    data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);

    NLODISFitter fitter(parameters);
    fitter.AddDataset(data);
    fitter.SetNLO(useNLO);
    fitter.SetSUB(useSUB);
    fitter.UseImprovedZ2Bound(useImprovedZ2Bound);
    fitter.UseConsistentlyBoundLoopTerm(useBoundLoop);
    fitter.SetCubaMethod(cubaMethod);

    cout << std::boolalpha;
    cout    << "=== Perturbative settings ===" << endl
            << "Settings: " << string_sub << " (scheme), " << string_bk << ", " << string_rc << endl
            << "Use ResumBK: " << !(config::LO_BK) << endl
            << "KinematicalConstraint / target eta0 BK: " << config::KINEMATICAL_CONSTRAINT << endl
            << "Running Coupling: (RC_LO):    " << config::RC_LO << " (0 fc, 1 parent, 4 balitsky, 6 guillaume)" << endl
            << "Running Coupling: (RESUM_RC): " << config::RESUM_RC << " (0 balitsky, 1 parent, 3 guillaume)" << endl
            << "Running Coupling: (RC_DIS):   " << nlodis_config::RC_DIS << " (0 fc, 1 parent, 2 guillaume)" << endl
            << "Use NLOimpact: " << useNLO << endl
            << "Use SUBscheme: " << useSUB << endl
            << "Use Sigma3: " << useSigma3 << endl
            << "Use improved Z2 bound: " << useImprovedZ2Bound << endl
            << "Use Z2 loop term: " << useBoundLoop << endl
            << "Cuba MC algorithm: " << cubaMethod << endl;
    cout << "=== Initial parameters ===" << endl;
    cout << parameters << endl;
    if(nlodis_config::VERBOSE) cout << "=== Starting fit ===" << endl;

    //MnMachinePrecision precc();
    //precc.SetPrecision(1e-2);
    // MnMinimize: use MIGRAD, if it fails, fall back to SIMPLEX
    // Optional 3rd argument, an unsigned int, set algorithm strategy: 0 = low (fast) , 1 = medium (def) , 2 = high (taxing)
    //MnMigrad fit(fitter, parameters, 0);
    MnMinimize fit(fitter, parameters);
    //MnSimplex fit(fitter,parameters);
    //MnScan fit(fitter, parameters);

    //fit.SetPrecision(0.01); // TODO Should this match BK solver acc?
    // minimize
    FunctionMinimum min = fit();
    // output
    std::cout<<"minimum: "<<min<<std::endl;

    cout << "=== Perturbative settings were ===" << endl;
    cout << std::boolalpha;
    cout    << "Use ResumBK: " << !(config::LO_BK) << endl
            << "KinematicalConstraint BK: " << config::KINEMATICAL_CONSTRAINT << endl
            << "Running Coupling: " << config::RC_LO << " (4 parent, 6 guillaume)" << endl
            << "Use NLOimpact: " << useNLO << endl
            << "Use SUBscheme: " << useSUB << endl
            << "Use Sigma3: " << useSigma3 << endl
            << "Use improved Z2 bound: " << useImprovedZ2Bound << endl
            << "Use Z2 loop term: " << useBoundLoop << endl
            << "Cuba MC algorithm: " << cubaMethod << endl;

    return 0;
}

//int errorscounter;
/*void ErrHandler(const char * reason,
                const char * file,
                int line,
                int gsl_errno)
{    // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    // 15: underflows

    if (gsl_errno == 15 or gsl_errno == 16) return;
    // Ugly hack, comes from the edges of the z integral in virtual_photon.cpp
    // Overflows come from IPsat::bint when it is done analytically
    // Hope is that these errors are handled correctly everywhere
    errorscounter++;
    std::cerr << file << ":"<< line <<": Error " << errorscounter << ": " <<reason
    << " (code " << gsl_errno << ")." << std::endl;
}
*/
