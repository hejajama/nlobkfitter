#include <cmath>
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

using namespace std;
using namespace ROOT::Minuit2;
//void ErrHandler(const char * reason,const char * file,int line,int gsl_errno);
int errors_counter;
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
    if (gsl_errno == 15) return; // underflow // safe to ignore?
    // if (gsl_errno == 16) return; // overflow
    // Ugly hack, comes from the edges of the z integral in virtual_photon.cpp
    // Overflows come from IPsat::bint when it is done analytically
    // Hope is that these errors are handled correctly everywhere
    errors_counter++;
    std::cerr << file << ":"<< line <<": Error " << errors_counter << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}

int main( int argc, char* argv[] )
{
    gsl_set_error_handler(&ErrHandlerCustom);

    // NLO DIS SIGMA_R COMPUTATION CONFIGS
    nlodis_config::USE_MASSES = false;

    nlodis_config::CUBA_EPSREL = 1e-3;
    // nlodis_config::CUBA_EPSREL = 5e-3; // highacc def1
    // nlodis_config::CUBA_MAXEVAL = 4e8;
    // nlodis_config::CUBA_MAXEVAL= 5e7;
    nlodis_config::CUBA_MAXEVAL= 2e7; // highacc def1
    // nlodis_config::CUBA_MAXEVAL= 1e10; // highacc big
    nlodis_config::MINR = 1e-6;
    // nlodis_config::MINR = 1e-7;
    // nlodis_config::MAXR = 12;
    nlodis_config::MAXR = 25;
    nlodis_config::PRINTDATA = true;
    bool useNLO = true;
    string cubaMethod = "vegas";
    // string cubaMethod = "suave";
    // string cubaMethod = "divonne";

    config::NO_K2 = true;  // Do not include numerically demanding full NLO part
    config::KINEMATICAL_CONSTRAINT = config::KC_NONE;

    config::VERBOSE = true;
    //config::RINTPOINTS = 512/4;
    //config::THETAINTPOINTS = 512/4;

    config::INTACCURACY = 10e-3;//0.02;
    // config::INTACCURACY = 2e-3;
    // config::INTACCURACY = 0.005; // used in the final fits
    //config::MINR = 1e-6;
    config::MINR = nlodis_config::MINR;
    //config::MAXR = 30;
    config::MAXR = nlodis_config::MAXR;
    //config::RPOINTS = 100;
    config::RPOINTS = 100;
    config::DE_SOLVER_STEP = 0.4; // Rungekutta step
    // config::DE_SOLVER_STEP = 0.8; // Rungekutta step

    // Constants
    config::NF=3;   // Only light quarks
    config::LAMBDAQCD = 0.241;


    MnUserParameters parameters;

    bool useSUB, useResumBK, useKCBK, useImprovedZ2Bound, useBoundLoop;
    bool calc_LO = true;
    bool useSigma3 = false;
    string helpstring = "Argument order: SCHEME BK RC useImprovedZ2Bound useBoundLoop Q C^2 X0 gamma Q0sq Y0 eta0\nsub/unsub/unsub+ resumbk/trbk/lobk parentrc/guillaumerc/fixedrc z2improved/z2simple z2boundloop/unboundloop";
    string string_sub, string_bk, string_rc;
    string ddis_calc;
    if (argc<2){ cout << helpstring << endl; return 0;}
    // Argv[0] is the name of the program

    ddis_calc = string(argv [1]);
    if (string(argv [1]) == "all"){
        useNLO = true;
    } else if (string(argv [1]) == "nlo"){
        useNLO = true;
        calc_LO = false;
    } else if (string(argv [1]) == "approx"){
        useNLO = false;
    } else {cout << helpstring << endl; return -1;}

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
    // nlodis_config::RC_SCALE_PRESC_DDIS = nlodis_config::DDIS_RC_SCALE_PRODUCT_COUPLING;
    // nlodis_config::RC_SCALE_PRESC_DDIS = nlodis_config::DDIS_RC_SCALE_SMALLER_DIS;
    nlodis_config::RC_SCALE_PRESC_DDIS = nlodis_config::DDIS_RC_SCALE_COMBINE_MEAN_DIS;
    if (string(argv [3]) == "parentrc" or string(argv [3]) == "pdrc"){
            config::RC_LO = config::PARENT_LO;
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
    } else {cout << helpstring << endl; return -1;}

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

    bool use_custom_prescription = false;
   string custom_presc = "none";
//    if (argc == 7){
//        use_custom_prescription = true;
//        custom_presc = argv[6];
//    }

    

    cout << std::boolalpha;
    cout    << "# === Perturbative settings ===" << endl
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
                << endl
            << "# config::INTACCURACY = " << config::INTACCURACY
                << ", config::RPOINTS = " << config::RPOINTS
                << ", config::DE_SOLVER_STEP = " << config::DE_SOLVER_STEP
                << ", config::{MINR, MAXR}, nlodis_config::{MINR, MAXR} = " << config::MINR << " " << config::MAXR << " " << nlodis_config::MINR << " " << nlodis_config::MAXR
                << endl;

    double qs0sqr;
    double alphas_scaling;
    double anomalous_dimension;
    double icx0_bk;
    double sigma02;

    if (argc == 6){
	    qs0sqr       = 0.2; //par[ parameters.Index("qs0sqr")];
        alphas_scaling     = 1.0; //par[ parameters.Index("alphascalingC2")]; // MATCH THIS IN IMPACTFACTOR ALPHA_S WHEN NLO
	    anomalous_dimension = 1.0; //par[ parameters.Index("anomalous_dimension")];
        icx0_bk = 0.01; //par[ parameters.Index("initialconditionX0")];
        sigma02 = 1.0;
	cout << "# NO INITIAL PARAMETERS PASSED, USING DEFAULTS" << endl;
    }

    if (argc >= 11){
    	qs0sqr		    = stod(argv[6]);
        alphas_scaling	    = stod(argv[7]);
	    anomalous_dimension = stod(argv[8]);
        icx0_bk		    = stod(argv[9]);
        sigma02		    = stod(argv[10]);

    }

    string dataname;
    if (argc >= 12){
        dataname = argv[11];
    }

    string string_ic = "cli_mode";
    /*
    *   PARAMS TO IC
    */
    double e_c          = 1.0; //par[ parameters.Index("e_c")];
    double icx0_nlo_impfac = 1.0; //par[ parameters.Index("initialconditionX0")];
    double initialconditionY0  = 0; //par[ parameters.Index("initialconditionY0")];
    double icTypicalPartonVirtualityQ0sqr  = 1.0; //par[ parameters.Index("icTypicalPartonVirtualityQ0sqr")];
    double qMass_light  = 0.14; // GeV --- doesn't improve fit at LO
    double qMass_u = 0.0023; // GeV, literature value
    double qMass_d = 0.0048; // GeV, literature value
    double qMass_s = 0.095; // GeV, literature value
    double qMass_charm = 1.35;
    double qMass_b = 4.180; // GeV, literature value
    bool useMasses = nlodis_config::USE_MASSES;
    bool useCharm = false;

    cout << "# === Initial parameters ===" << endl;
    cout << "# "
         << "qs0sqr=" << qs0sqr
         << ", alphas_scaling=" << alphas_scaling
         << ", gamma=" << anomalous_dimension
         << ", e_c=" << e_c
         << ", icx0_nlo_impfac=" << icx0_nlo_impfac
         << ", icx0_bk=" << icx0_bk
         << ", icY0=" << initialconditionY0
         << ", icTypPartonVirtQ0sqr=" << icTypicalPartonVirtualityQ0sqr
         << ", sigma02=" << sigma02
         << endl;

    /*
    // ***Solve BK***
    */
    double eta0 = 0;

    double min_xbj = 1e-7;
    // NLODISSIGMAR IC + SOLVER CODE
    MV ic;                                            // Initial condition
    ic.SetQsqr(qs0sqr);
    ic.SetAnomalousDimension(anomalous_dimension);
    ic.SetE(e_c);                                     // e_c of MVe parametrization

    Dipole dipole(&ic);
    dipole.SetX0(icx0_bk);
    BKSolver solver(&dipole);
    double maxy = std::log(icx0_bk/(min_xbj)) + initialconditionY0; // divisor=smallest HERA xbj in Q^2 range (1E-05)?
    
    double maxq2=100;     // Todo, take from actual data!
    if (useImprovedZ2Bound)
        maxy += std::log(maxq2 / icTypicalPartonVirtualityQ0sqr);

    // cout << "=== Solving BK ===" << endl;

    solver.SetAlphasScaling(alphas_scaling);
    // solver.SetEta0(par[ parameters.Index("eta0")]);
    solver.SetX0(icx0_bk);
    solver.SetICX0_nlo_impfac(icx0_nlo_impfac);
    solver.SetICTypicalPartonVirtualityQ0sqr(icTypicalPartonVirtualityQ0sqr);



    AmplitudeLib* DipoleAmplitude_ptr; // Forward declaration of the dipole object to be initialized from a file or solved data.
    string dipole_basename = "./out/dipoles/dipole";
    string dipole_filename = dipole_basename
                             + "_" + string_bk
                             + "_" + string_rc
                             + "_x0bk" + std::to_string(icx0_bk)
                             + "_qs0sqr" + std::to_string(qs0sqr)
                             + "_asC^2" + std::to_string(alphas_scaling)
                             + "_gamma" + std::to_string(anomalous_dimension)
                             + "_ec" + std::to_string(e_c)
                             + "_eta0" + std::to_string(eta0)
                             + "_maxy" + std::to_string(maxy)
                             + "_euler" + std::to_string(config::EULER_METHOD)
                             + "_step" + std::to_string(config::DE_SOLVER_STEP)
                             + "_rpoints" + std::to_string(config::RPOINTS)
                             + "_rminmax" + std::to_string(config::MINR) + "--" + std::to_string(config::MAXR)
                             + "_intacc" + std::to_string(config::INTACCURACY) ;
    
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
    } else {
        solver.Solve(maxy);     // Solve up to maxy since specified dipole datafile was not found.
        solver.GetDipole()->Save(dipole_filename);
        cout << "# Saved dipole to file: "<< dipole_filename << endl;
        DipoleAmplitude_ptr = new AmplitudeLib(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());
    }   
    // solver.GetDipole()->Save("./out/dipoles/dipole_lobk_fc_RK_step0.2_rpoints100_rmin1e-6rmax50_INTACC10e-3.dat");
    // cout << "Saved dipole to file, exiting." << endl;
    // exit(0);

    // Just solve dipole everytime:
    // solver.Solve(maxy);     // Solve up to maxy since specified dipole datafile was not found.
    // solver.GetDipole()->Save(dipole_filename);
    // cout << "# Saved dipole to file: "<< dipole_filename << endl;
    // DipoleAmplitude_ptr = new AmplitudeLib(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());

    // Give solution to the AmplitudeLib object
    // AmplitudeLib DipoleAmplitude(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());
    // AmplitudeLib DipoleAmplitude("./data/paper1dipole/pap1_fcBK_MV.dat"); // pap1_fcBK_MV.dat, pap1_rcBK_MV_parent.dat
    AmplitudeLib DipoleAmplitude(*DipoleAmplitude_ptr);
    // delete DipoleAmplitude_ptr;
    DipoleAmplitude.SetInterpolationMethod(LINEAR_LINEAR);
    DipoleAmplitude.SetX0(icx0_bk);
    DipoleAmplitude.SetOutOfRangeErrors(false);
    AmplitudeLib *DipolePointer = &DipoleAmplitude;

    ComputeSigmaR SigmaComputer(DipolePointer);
    SigmaComputer.SetX0(icx0_nlo_impfac);
    SigmaComputer.SetX0_BK(icx0_bk);
    SigmaComputer.SetY0(initialconditionY0);
    SigmaComputer.SetQ0Sqr(icTypicalPartonVirtualityQ0sqr);
    SigmaComputer.SetQuarkMassLight(qMass_light);
    SigmaComputer.SetQuarkMassCharm(qMass_charm);
    SigmaComputer.SetAlphasScalingC2(alphas_scaling);

    // Set running coupling and rapidity function pointters
    SigmaComputer.MetaPrescriptionSetter();
    
    // CUBA Monte Carlo integration library algorithm setter
    SigmaComputer.SetCubaMethod(cubaMethod);

    cout << "# === Computing Reduced Cross sections ===" << endl;

            // Print column titles.
            if(!useMasses){
            #pragma omp critical
            cout    << setw(15) << "# xpom"        << " "
                    << setw(15) << "Q^2"          << " "
                    << setw(15) << "beta"         << " "
                    << setw(15) << "xpom*FL_qq"        << " "
                    << setw(15) << "xpom*FL_qqg"       << " "
                    << setw(15) << "xpom*FT_qq"        << " "
                    << setw(15) << "xpom*FT_qqg"       << " "
                    << setw(15) << "xpom*FT_M^2>>"     << " "
                    << setw(15) << "xpom*FT_Q^2>>"     << " "
                    << endl;
                    }


    double icQ = 1.0;
    //double icx0 = icx0_bk;
    double icx0 = 0.01;

    std::vector< std::tuple<int, int, int> > coordinates;

    // for (int k=0; k<1; k+=1) // one BETA bin
    for (int k=0; k<10; k+=1) // BETA over a range
    {
        // for (int i=0; i<=20; i+=17)  // Q^2 = {1,50}
        // for (int i=0; i<=20; i+=1)  // Q^2 in [1,100]
        // for (int i=0; i<=20; i+=10)  // Q^2 in [1,100]
        // for (int i=1; i<=17; i+=4)  // Q^2 in [1,100]
        // for (int i=10; i<=20; i+=11)  // Q^2 = 10 (5)
        for (int i=20; i<=20; i+=11)  // Q^2 = 100
        // for (int i=0; i<=1; i++)
        {
            // for (int j=0; j<=17; j++)  // xbj in [5.62341e-07, 1e-2]
            //for (int j=4; j<=12; j+=8)  // xbj = {1e-3, 1e-5}
            // for (int j=1; j<=17; j+=8)  // xbj = {~1e-2, ~1e-4, ~1e-6} // LHEC predictions for x0bk=0.01
            for (int j=8; j<=10; j+=8)  // xbj = {1e-3}
            // for (int j=1; j<=9; j+=8)      // xbj = {~1e-2, ~1e-4} // LHEC predictions for x0bk=0.01
            // for (int j=2; j<=8; j+=2)  // xbj = {1e-6} // LHEC predictions for x0bk=0.01
            // for (int j=0; j<=1; j++)
            {
                if (!((i == 0 or i == 17) or (j == 1 or j == 4 or j == 8 or j == 12 or j == 17))) { continue; }
                coordinates.emplace_back(i,j,k);
            }
        }
    }

    #pragma omp parallel for
    for (size_t k=0; k< coordinates.size(); k++)
        {
        int i,j,l;
        std::tie(i,j,l) = coordinates[k];

        if (j==0 and cubaMethod=="suave"){continue;}
        // double Q = 1.0*pow(10,(double)i/20.0);
        double Q = 1.0*pow(5,(double)i/20.0);
        double xpom = icx0/pow(10,(double)j/4.0);
        // double beta = 0.95/pow(10,(double)l/4.0);
        double beta = 0.05 + (double)l * 0.1;
        // double beta = 0.5;
        
        double FL_qq=0, FL_qqg=0;
        double FT_qq=0, FT_qqg=0, FT_qqg_wusthof=0, FT_qqg_muniershoshi=0;

        int calccount=0;
        if (!useNLO && !useMasses) // Compute reduced cross section using leading order impact factors
        {
            FL_qq = SigmaComputer.diff_lo_xpom_FL(Q,xpom,beta);
            // cout << FL_qq << endl;
            FT_qq = SigmaComputer.diff_lo_xpom_FT(Q,xpom,beta);
            // cout << FT_qq << endl;
            FT_qqg_muniershoshi = SigmaComputer.diff_nlo_xpom_FT_qqbarg_largeM(Q,xpom,beta);
            // cout << FT_qqg_muniershoshi << endl;
            FT_qqg_wusthof = SigmaComputer.diff_nlo_xpom_FT_qqbarg_largeQsq(Q,xpom,beta);
            // cout << FT_qqg_wusthof << endl;
            ++calccount;
        }

        if (useNLO) // UNSUB SCHEME Full NLO impact factors for reduced cross section
        {
            if (calc_LO){
                FL_qq = SigmaComputer.diff_lo_xpom_FL(Q,xpom,beta);
                // cout << FL_qq << endl;
                FT_qq = SigmaComputer.diff_lo_xpom_FT(Q,xpom,beta);
                // cout << FT_qq << endl;
                FL_qqg = SigmaComputer.diff_nlo_xpom_FL_qqbarg(Q,xpom,beta);
                // cout << FL_qqg << endl;
                FT_qqg = SigmaComputer.diff_nlo_xpom_FT_qqbarg(Q,xpom,beta);
                // cout << FT_qqg << endl;
                FT_qqg_muniershoshi = SigmaComputer.diff_nlo_xpom_FT_qqbarg_largeM(Q,xpom,beta);
                // cout << FT_qqg_muniershoshi << endl;
                FT_qqg_wusthof = SigmaComputer.diff_nlo_xpom_FT_qqbarg_largeQsq(Q,xpom,beta);
                // cout << FT_qqg_wusthof << endl;
                ++calccount;}
            else{
                FL_qqg = SigmaComputer.diff_nlo_xpom_FL_qqbarg(Q,xpom,beta);
                FT_qqg = SigmaComputer.diff_nlo_xpom_FT_qqbarg(Q,xpom,beta);
            }
            if (useMasses){
                ++calccount;}
        }

        if (calccount>1)
        {
            cerr << "ERROR: Multiple computations. abort." << "count="<< calccount << endl;
            exit(1);
        }

        // Output for plotting
        bool print_ratios = false;
        if(!useMasses){
        #pragma omp critical
        cout    << setw(15) << xpom          << " "
                << setw(15) << Q*Q           << " "
                << setw(15) << beta          << " "
                << setw(15) << sigma02*FL_qq         << " "
                << setw(15) << sigma02*FL_qqg        << " "
                << setw(15) << sigma02*FT_qq         << " "
                << setw(15) << sigma02*FT_qqg        << " "
                << setw(15) << sigma02*FT_qqg_muniershoshi       << " "
                << setw(15) << sigma02*FT_qqg_wusthof            << " "
                << endl;
                }
    }
    return 0;
}
