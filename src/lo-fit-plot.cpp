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

    string method_string = "suave";
    double cubaeps = 1e-4;
    double cubamaxev = 2e7;
    double minuit_prec = 1e-10;
    // if (argc >= 19){
	// method_string = string(argv[15]);
	// cubaeps = stod( argv[16] );
	// cubamaxev = stod( argv[17] );
	// minuit_prec = stod( argv[18] );
    // }
    bool fit_qs0sqr = false, fit_csqr = false, fit_gamma = false, fit_ec = false;
    // if (argc == 20){
    //     string fit_par_string = argv[19];
    //     if (fit_par_string.find('0') != std::string::npos){
    //         fit_qs0sqr = true;
    //     }
    //     if (fit_par_string.find('1') != std::string::npos){
    //         fit_csqr = true;
    //     }
    //     if (fit_par_string.find('2') != std::string::npos){
    //         fit_gamma = true;
    //     }
    //     if (fit_par_string.find('3') != std::string::npos){
    //         fit_ec = true;
    //     }
    // }

    // NLO DIS SIGMA_R COMPUTATION CONFIGS
    nlodis_config::CUBA_EPSREL = cubaeps;
    nlodis_config::CUBA_MAXEVAL= cubamaxev;
    nlodis_config::MINR = 1e-6;
    nlodis_config::MAXR = 50;
    nlodis_config::PRINTDATA = true;
    bool useNLO = false;
    cout << "Running LO fit plot tool, useNLO is: " << useNLO << endl;
    string cubaMethod = method_string;

    // Constants
    config::NF=3;   // Only light quarks
    config::LAMBDAQCD = 0.241;
    
    config::VERBOSE = true;
    config::KSUB = 0.65;  // Optimal value for K_sub
    config::NO_K2 = true;  // Do not include numerically demanding full NLO part
    config::KINEMATICAL_CONSTRAINT = config::KC_NONE;

    config::INTACCURACY = 5e-3;//0.02;
    // config::INTACCURACY = 1e-3;//0.02;
    config::MINR = 1e-6;
    config::MAXR = 50;
    config::RPOINTS = 200;
    // config::DE_SOLVER_STEP = 0.4; // Rungekutta step
    config::DE_SOLVER_STEP = 0.1; // Rungekutta step
    // config::DE_SOLVER_STEP = 0.05; // Euler step


    // READING RUN CONFIGURATION FROM THE STDIN
    // bool useSUB, useResumBK, useKCBK, useImprovedZ2Bound, useBoundLoop, useSigma3;
    bool useSUB = false;
    bool useResumBK = false;
    bool useImprovedZ2Bound = false;
    bool useBoundLoop = false;
    bool useSigma3 = false;
    // string helpstring = "Argument order: SCHEME BK RC useImprovedZ2Bound useBoundLoop [Qs0 C^2 gamma] X0_if X0_bk e_c Q0sq Y0 eta0 CUBA_MTHD CUBA_EPS CUBA_MEVAL MINUIT_PREC\nsub/unsub/unsub+ resumbk/trbk/lobk parentrc/guillaumerc/fixedrc z2improved/z2simple z2boundloop/unboundloop";
    string helpstring = "Configure in lo-fit-plot.cpp";
    // string string_sub, string_bk, string_rc;
    string string_sub = "N/A";
    // if (argc<2){ cout << helpstring << endl; return 0;}
    // Argv[0] is the name of the program

    string string_bk = "lobk";
    string string_rc = "sdrc";
    // string string_ic_preset = "hmtl-mv";
    // string string_ic_preset = "hmtl-mvgamma";
    string string_ic_preset = "hmtl-mvec";
    string fit_par_string = ""; // empty string doesn't fit, 0 = qs0sqr, 1 = C^2, 2 = gamma, 3 = ec


    // string_bk = string(argv [2]);
    if (string(string_bk) == "resumbk"){
            config::EULER_METHOD = false;    // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = true;       // Resum doulbe logs
            config::RESUM_SINGLE_LOG = true; // Resum single logs
            config::KSUB = 0.65;             // Optimal value for K_sub
            config::NO_K2 = true;            // Do not include numerically demanding full NLO part
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_RESUM;
    }else if (string(string_bk) == "kcbk"){
            config::EULER_METHOD = true;     // Kinematical constraint requires this
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_BEUF_K_PLUS;
            config::DE_SOLVER_STEP = 0.05;  //0.02; // Euler method requires smaller step than RungeKutta!
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_KCBK_BEUF;
    }else if (string(string_bk) == "trbk"){  // Target Rapidity BK
            config::EULER_METHOD = true;     // Kinematical constraint requires this
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_EDMOND_K_MINUS;
            config::DE_SOLVER_STEP = 0.05;  //0.02; // Euler method requires smaller step than RungeKutta!
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_TRBK_EDMOND;
            nlodis_config::TRBK_RHO_PRESC = nlodis_config::TRBK_RHO_QQ0;
    }else if (string(string_bk) == "lobk"){
            config::EULER_METHOD = false;   // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_NONE;
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_LOBK_EXPLICIT;
    } else {cout << helpstring << endl; return -1;}

    // string_rc = string(argv [3]);
    if (string(string_rc) == "parentrc" or string(string_rc) == "pdrc"){
            config::RC_LO = config::PARENT_LO;
            config::RESUM_RC = config::RESUM_RC_PARENT;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_PARENT;
    } else if (string(string_rc) == "guillaumerc" or string(string_rc) == "gbrc"){
            config::RC_LO = config::GUILLAUME_LO;
            config::RESUM_RC = config::RESUM_RC_GUILLAUME;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_GUILLAUME;
    } else if (string(string_rc) == "fixedrc" or string(string_rc) == "fc"){
            config::RC_LO = config::FIXED_LO;
            config::RESUM_RC = config::RESUM_RC_FIXED;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_FIXED;
    } else if (string(string_rc) == "smallestrc" or string(string_rc) == "sdrc"){
            config::RC_LO = config::BALITSKY_LO;
            config::RESUM_RC = config::RESUM_RC_SMALLEST;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_SMALLEST;
    } else {cout << helpstring << endl; return -1;}


    double icqs0sq, iccsq, icx0_if, icx0_bk, ic_ec, icgamma, icQ0sq, icY0, icEta0;
    // reading fit initial condition parameters:
    if (string_ic_preset == "hmtl-mv") {
        // use LO fit as IC
        icqs0sq   = 0.104;
        iccsq     = 14.5;
        icgamma   = 1.0;
        ic_ec     = 1.0;

        icx0_if   = 1.0;
        icx0_bk   = 0.01;
        icQ0sq    = 1.0;
        icY0      = 0.0;
        icEta0    = 0.0;
    } else if (string_ic_preset == "hmtl-mvgamma") {
        // use LO fit as IC
        icqs0sq   = 0.165;
        iccsq     = 6.35;
        icgamma   = 1.135;
        ic_ec     = 1.0;

        icx0_if   = 0.01;
        icx0_bk   = 0.01;
        icQ0sq    = 1.0;
        icY0      = 0.0;
        icEta0    = 0.0;
    } else if (string_ic_preset == "hmtl-mvec") {
        // use LO fit as IC
        icqs0sq   = 0.060;
        iccsq     = 7.2;
        icgamma   = 1.0;
        ic_ec     = 18.9;

        icx0_if   = 1.0;
        icx0_bk   = 0.01;
        icQ0sq    = 1.0;
        icY0      = 0.0;
        icEta0    = 0.0;
    } else {
        cout << "Bad fit ic preset!" << endl;
        exit(1);
    }

    if (fit_par_string.find('0') != std::string::npos){
        fit_qs0sqr = true;
    }
    if (fit_par_string.find('1') != std::string::npos){
        fit_csqr = true;
    }
    if (fit_par_string.find('2') != std::string::npos){
        fit_gamma = true;
    }
    if (fit_par_string.find('3') != std::string::npos){
        fit_ec = true;
    }

    MnUserParameters parameters;
    // Fit parameters, first value is starting value, second is uncertainty
        if (fit_qs0sqr){
            parameters.Add("qs0sqr",                icqs0sq, 0.01);
        }else{
            parameters.Add("qs0sqr",                icqs0sq);
        }

        if (fit_csqr){
            parameters.Add("alphascalingC2",        iccsq, 5.0);
        }else{
            parameters.Add("alphascalingC2",        iccsq);
        }
        
        if (fit_gamma){
            parameters.Add("anomalous_dimension",   icgamma, 0.2);
        }else{
            parameters.Add("anomalous_dimension",   icgamma);
        }
                
        parameters.Add("icx0_nlo_impfac",       icx0_if );
        parameters.Add("icx0_bk",               icx0_bk );
        if (fit_ec){
            parameters.Add("e_c",                   ic_ec , 0.1);
        }else{
            parameters.Add("e_c",                   ic_ec );
        }        
        
        parameters.Add("icTypicalPartonVirtualityQ0sqr", icQ0sq );
        parameters.Add("initialconditionY0",    icY0 );
        parameters.Add("eta0", icEta0 );
    //
    // Set limits
    //
        // parameters.SetLimits("qs0sqr",              0.01, 0.2);
        // parameters.SetLimits("alphascalingC2",      0.1, 200.0);
        // parameters.SetLimits("anomalous_dimension", 0.9 , 3.0);
        //parameters.SetLimits("icx0_nlo_impfac",	    0.01 , 10.0);



    Data data;
    // data.SetMinQsqr(0.75);
    data.SetMinQsqr(1.0);
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
    cout    << "# === Perturbative settings ===" << endl
            << "# Settings: " << string_sub << " (scheme), " << string_bk << ", " << string_rc << endl
            << "# Use LOBK (DL,SL==false): " << (!(config::RESUM_DLOG) 
                                    and !(config::RESUM_SINGLE_LOG)) << endl
            << "# Use ResumBK (DL,SL==true,KC_NONE): " << ((config::RESUM_DLOG) 
                                    and (config::RESUM_SINGLE_LOG)
                                    and (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)) << endl
            << "# KinematicalConstraint / target eta0 BK: " << config::KINEMATICAL_CONSTRAINT << " (0 BEUF_K_PLUS, 1 EDMOND_K_MINUS, 2 NONE)" << endl
            << "# Running Coupling: (RC_LO):    " << config::RC_LO << " (0 fc, 1 parent, 4 balitsky, 6 guillaume)" << endl
            << "# Running Coupling: (RESUM_RC): " << config::RESUM_RC << " (0 fc, 1 balitsky, 2 parent, 4 guillaume)" << endl
            << "# Running Coupling: (RC_DIS):   " << nlodis_config::RC_DIS << " (0 fc, 1 parent, 2 guillaume)" << endl
            << "# Use NLOimpact: " << useNLO << endl
            << "# Use SUBscheme: " << useSUB << endl
            << "# Use Sigma3: " << useSigma3 << endl
            << "# Use improved Z2 bound: " << useImprovedZ2Bound << endl
            << "# Use Z2 loop term: " << useBoundLoop << endl
            << "# Cuba MC: " << cubaMethod
                << ", Cuba eps = " << nlodis_config::CUBA_EPSREL
                << ", Cuba maxeval = " << (float)nlodis_config::CUBA_MAXEVAL
                << endl
	    << "# Minuit Precision = " << minuit_prec
		<< endl
            << "# config::INTACCURACY = " << config::INTACCURACY
                << ", config::RPOINTS = " << config::RPOINTS
                << ", config::DE_SOLVER_STEP = " << config::DE_SOLVER_STEP
                << ", config::{MINR, MAXR}, nlodis_config::{MINR, MAXR} = " << config::MINR << " " << config::MAXR << " " << nlodis_config::MINR << " " << nlodis_config::MAXR
                << endl;
    cout << "=== Initial parameters ===" << endl;
    cout << parameters << endl;
    if(nlodis_config::VERBOSE) cout << "=== Starting fit ===" << endl;

    // MnMinimize: use MIGRAD, if it fails, fall back to SIMPLEX
    // Optional 3rd argument, an unsigned int, set algorithm strategy: 0 = low (fast) , 1 = medium (def) , 2 = high (taxing)
    //MnMigrad fit(fitter, parameters, 0);
    MnMinimize fit(fitter, parameters);
    //MnSimplex fit(fitter,parameters);
    //MnScan fit(fitter, parameters);

    if (minuit_prec != 0){
        fit.SetPrecision(minuit_prec);
    }
    // minimize
    FunctionMinimum min = fit();
    // output
    std::cout<<"minimum: "<<min<<std::endl;

    cout << "=== Perturbative settings were ===" << endl;
    cout << std::boolalpha;
    cout    << "# === Perturbative settings ===" << endl
            << "# Settings: " << string_sub << " (scheme), " << string_bk << ", " << string_rc << endl
            << "# Use LOBK (DL,SL==false): " << (!(config::RESUM_DLOG) 
                                    and !(config::RESUM_SINGLE_LOG)) << endl
            << "# Use ResumBK (DL,SL==true,KC_NONE): " << ((config::RESUM_DLOG) 
                                    and (config::RESUM_SINGLE_LOG)
                                    and (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)) << endl
            << "# KinematicalConstraint / target eta0 BK: " << config::KINEMATICAL_CONSTRAINT << " (0 BEUF_K_PLUS, 1 EDMOND_K_MINUS, 2 NONE)" << endl
            << "# Running Coupling: (RC_LO):    " << config::RC_LO << " (0 fc, 1 parent, 4 balitsky, 6 guillaume)" << endl
            << "# Running Coupling: (RESUM_RC): " << config::RESUM_RC << " (0 fc, 1 balitsky, 2 parent, 4 guillaume)" << endl
            << "# Running Coupling: (RC_DIS):   " << nlodis_config::RC_DIS << " (0 fc, 1 parent, 2 guillaume)" << endl
            << "# Use NLOimpact: " << useNLO << endl
            << "# Use SUBscheme: " << useSUB << endl
            << "# Use Sigma3: " << useSigma3 << endl
            << "# Use improved Z2 bound: " << useImprovedZ2Bound << endl
            << "# Use Z2 loop term: " << useBoundLoop << endl
            << "# Cuba MC: " << cubaMethod
                << ", Cuba eps = " << nlodis_config::CUBA_EPSREL
                << ", Cuba maxeval = " << (float)nlodis_config::CUBA_MAXEVAL
                << endl
	    << "# Minuit Precision = " << minuit_prec
		<< endl
            << "# config::INTACCURACY = " << config::INTACCURACY
                << ", config::RPOINTS = " << config::RPOINTS
                << ", config::DE_SOLVER_STEP = " << config::DE_SOLVER_STEP
                << ", config::{MINR, MAXR}, nlodis_config::{MINR, MAXR} = " << config::MINR << " " << config::MAXR << " " << nlodis_config::MINR << " " << nlodis_config::MAXR
                << endl;

    return 0;
}
