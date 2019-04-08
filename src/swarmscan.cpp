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
#include "nlodissigmar.hpp"
#include "nlodis_config.hpp"

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
//void ErrHandler(const char * reason,const char * file,int line,int gsl_errno);
int errors_mmyiss;
void ErrHandlerCustom(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{   // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    // 15: underflows

    if (gsl_errno == 15 or gsl_errno == 16) return;
    // Ugly hack, comes from the edges of the z integral in virtual_photon.cpp
    // Overflows come from IPsat::bint when it is done analytically
    // Hope is that these errors are handled correctly everywhere
    errors_mmyiss++;
    std::cerr << file << ":"<< line <<": Error " << errors_mmyiss << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}

int main( int argc, char* argv[] )
{
    gsl_set_error_handler(&ErrHandlerCustom);

        // NLO DIS SIGMA_R COMPUTATION CONFIGS
        nlodis_config::CUBA_EPSREL = 2e-2;
        nlodis_config::CUBA_MAXEVAL= 1e8;
        nlodis_config::PRINTDATA = true;
        bool useNLO = true;
        //bool useSUB = true;               // Set by a command line switch in swarmscan
        //bool useImprovedZ2Bound = true;   // Set by a command line switch in swarmscan
        //bool useBoundLoop = true;         // Set by a command line switch in swarmscan
        string cubaMethod = "suave";
        //bool useResumBK = true;
        //bool useKCBK = false;
        //if (useResumBK == true and useKCBK == true) {cout << "Both ResumBK and KCBK enabled, exitting." << endl; return -1;}

        //config::RC_LO = config::FIXED_LO;
        config::RC_LO = config::BALITSKY_LO; // Balitsky running coupling for LO kernel
        //config::RC_LO = config::GUILLAUME_LO;
        config::RESUM_RC = config::RESUM_RC_PARENT; // Parent dipole in the resummation
        config::RESUM_DLOG = true; // Resum doulbe logs
        config::RESUM_SINGLE_LOG = true; // Resum single logs
        config::LO_BK = true;  // Solve LO BK with running coupling, overrides RESUM settings
        config::KSUB = 0.65;  // Optimal value for K_sub
        config::NO_K2 = true;  // Do not include numerically demanding full NLO part
        config::INTACCURACY = 1e-2;//0.02;
        config::MINR = 1e-5;
        config::MAXR = 30;
        config::RPOINTS = 100;
        config::DE_SOLVER_STEP = 0.4; // Rungekutta step

        // Constants
        config::NF=3;   // Only light quarks
        config::LAMBDAQCD = 0.241;


    Data data;
    data.SetMinQsqr(0.75);
    data.SetMaxQsqr(50);
    data.SetMaxX(0.01);
    data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);

    MnUserParameters parameters;

    bool useSUB, useResumBK, useKCBK, useImprovedZ2Bound, useBoundLoop, useSigma3;
    string helpstring = "Argument order: SCHEME BK RC useImprovedZ2Bound useBoundLoop Q C^2 X0 gamma Q0sq Y0";
    if (argc<2){ cout << helpstring << endl; return 0;}
    // Argv[0] is the name of the program

    if (string(argv [1]) == "sub"){
      useSUB = true;
    } else if (string(argv [1]) == "unsub"){
      useSUB = false;
      useSigma3 = false;
    } else if (string(argv [1]) == "unsub+"){
      useSUB = false;
      useSigma3 = true;
    } else {cout << helpstring << endl; return -1;}

    if (string(argv [2]) == "resumbk"){
            config::LO_BK = false;  // Solve LO BK with running coupling, overrides RESUM settings
    } else if (string(argv [2]) == "kcbk"){
            config::LO_BK = true;  // Solve LO BK with running coupling, overrides RESUM settings
            config::EULER_METHOD            = true;        // Kinematical constraint requires this
            config::KINEMATICAL_CONSTRAINT  = true;
            config::DE_SOLVER_STEP = 0.08; //0.02; // Euler method requires smaller step than RungeKutta!
    }else if (string(argv [2]) == "lobk"){
            config::LO_BK = true;
    } else {cout << helpstring << endl; return -1;}

    if (string(argv [3]) == "parentrc"){
            config::RC_LO = config::BALITSKY_LO; 
    } else if (string(argv [3]) == "guillaumerc"){
            config::RC_LO = config::GUILLAUME_LO;
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
    double QAdditN   = stod( argv [argi] ); argi++;
    double CPowerN   = stod( argv [argi] ); argi++;
    double X0calc    = stod( argv [argi] ); argi++;
    double gammacalc = stod( argv [argi] ); argi++;
    double icQ0sq    = stod( argv [argi] ); argi++;
    double icY0      = stod( argv [argi] );
    icY0 *= std::log(10.0/X0calc);


    // Constructing Q and C^2 grids
    double Qdown = 0;
    double Qup = 0.2;
    double QStep = (Qup-Qdown)/200.;
    double Qcalc = Qdown + QAdditN * QStep;

    double Cdown = 0.1;
    double Cup = 5.0;
    double CStep = pow(Cup/Cdown , CPowerN/240.);
    double Ccalc = Cdown * CStep;


    cout << Qcalc << " " << Ccalc << endl;

	parameters.Add("qs0sqr",		        Qcalc);
        parameters.Add("fitsigma0",		        1.0); // 1mb = 2.568 GeV²
        parameters.Add("alphascalingC2",	    Ccalc);
        parameters.Add("e_c",                   1.0 );
        parameters.Add("anomalous_dimension",   gammacalc );
        parameters.Add("initialconditionX0",    X0calc ); // fixed initialx0 computations were done with 0.01
        parameters.Add("initialconditionY0",    icY0 );
        parameters.Add("icTypicalPartonVirtualityQ0sqr", icQ0sq );

    NLODISFitter fitter(parameters);
    fitter.AddDataset(data);
    fitter.SetNLO(useNLO);
    fitter.SetSUB(useSUB);
    fitter.SetSigma3(useSigma3);
    fitter.UseImprovedZ2Bound(useImprovedZ2Bound);
    fitter.UseConsistentlyBoundLoopTerm(useBoundLoop);
    fitter.SetCubaMethod(cubaMethod);

    cout << "=== Perturbative settings ===" << endl;
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
    cout << "=== Initial parameters ===" << endl;
    cout << parameters << endl;
    if(nlodis_config::VERBOSE) cout << "=== Starting fit ===" << endl;

    MnMinimize fit(fitter, parameters);
    FunctionMinimum min = fit();
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
