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
    if (gsl_errno == 15 or gsl_errno == 16) return;
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
    nlodis_config::USE_MASSES = true;


    nlodis_config::CUBA_EPSREL = 1e-3;
    //nlodis_config::CUBA_EPSREL = 5e-3; // highacc def1
    nlodis_config::CUBA_MAXEVAL = 6e7;
    // nlodis_config::CUBA_MAXEVAL= 5e7; // highacc def1
    //nlodis_config::MINR = 1e-6;
    nlodis_config::MINR = 1e-6;
    nlodis_config::MAXR = 30;
    nlodis_config::PRINTDATA = true;
    bool useNLO = true;
    bool computeNLO = useNLO;
    string cubaMethod = "vegas";
    // string cubaMethod = "suave";

    config::NO_K2 = true;  // Do not include numerically demanding full NLO part
    config::KINEMATICAL_CONSTRAINT = config::KC_NONE;

    config::VERBOSE = true;
    //config::RINTPOINTS = 512/4;
    //config::THETAINTPOINTS = 512/4;

    // config::INTACCURACY = 10e-3;//0.02;
    config::INTACCURACY = 10e-3;
    //config::MINR = 1e-6;
    config::MINR = 1e-6;
    //config::MAXR = 30;
    config::MAXR = 30;
    //config::RPOINTS = 100;
    config::RPOINTS = 100;
    config::DE_SOLVER_STEP = 0.4; // Rungekutta step
    // config::DE_SOLVER_STEP = 0.8; // Rungekutta step

    // Constants
    config::NF=3;   // Only light quarks
    config::LAMBDAQCD = 0.241;


    Data data;
    data.SetMinQsqr(0.75);
    data.SetMaxQsqr(100);
    data.SetMaxX(0.01);
    data.LoadData("./data/h1_fl_averageddata_1312.4821tab6.txt", FL);
    vector<Data*> datasets;
    datasets.push_back(&data);

    MnUserParameters parameters;

    bool useSUB, useResumBK, useKCBK, useImprovedZ2Bound, useBoundLoop;
    bool useSigma3 = false;
    string helpstring = "Argument order: SCHEME BK RC useImprovedZ2Bound useBoundLoop Q C^2 X0 gamma Q0sq Y0 eta0\nsub/unsub/unsub+ resumbk/trbk/lobk parentrc/guillaumerc/fixedrc z2improved/z2simple z2boundloop/unboundloop";
    string string_sub, string_bk, string_rc;
    if (argc<2){ cout << helpstring << endl; return 0;}
    // Argv[0] is the name of the program

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
            config::RC_LO = config::BALITSKY_LO;
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

    double qs0sqr;
    double alphas_scaling;
    double anomalous_dimension;
    double icx0_bk;
    double sigma02;
    double qmass;

    if (argc == 6){
	    qs0sqr       = 0.2; //par[ parameters.Index("qs0sqr")];
        alphas_scaling     = 1.0; //par[ parameters.Index("alphascalingC2")]; // MATCH THIS IN IMPACTFACTOR ALPHA_S WHEN NLO
	    anomalous_dimension = 1.0; //par[ parameters.Index("anomalous_dimension")];
        icx0_bk = 0.01; //par[ parameters.Index("initialconditionX0")];
        sigma02 = 1.0;
	cout << "# NO INITIAL PARAMETERS PASSED, USING DEFAULTS" << endl;
    }

    if (argc == 11){
    	qs0sqr		    = stod(argv[6]);
        alphas_scaling	    = stod(argv[7]);
	    anomalous_dimension = stod(argv[8]);
        icx0_bk		    = stod(argv[9]);
        sigma02		    = stod(argv[10]);
    }
    double icqs0sq, iccsq, icx0_if, ic_ec=1.0, icgamma, icQ0sq, icY0, icEta0; 
    double Cdown, Cup, CStep;
    if (argc >= 15){
        icqs0sq   = stod( argv [6] );
        iccsq     = stod( argv [7] );
        anomalous_dimension   = stod( argv [8] );
        icx0_if   = stod( argv [9] );
        icx0_bk   = stod( argv [10] );
        ic_ec     = stod( argv [11] );
        icQ0sq    = stod( argv [12] );
        icY0      = stod( argv [13] );
        icEta0    = stod( argv [14] );
        sigma02 = stod( argv [15] );
        qmass = stod( argv[16] );

        // Constructing Q and C^2 grids
        qs0sqr = icqs0sq/1000.; // input in milli-GeV² integer strings

        Cdown = 0.1;
        Cup = 5.0;
        CStep = pow(Cup/Cdown , iccsq/240.);
        alphas_scaling = Cdown * CStep; // compute C² on a logarithmic grid
    }

    string string_ic = "cli_mode";
    /*
    *   PARAMS TO IC
    */
    if (string_ic == "default"){
    cout << "this should not trigger." << endl;
    qs0sqr       = 0.2; //par[ parameters.Index("qs0sqr")];
    alphas_scaling     = 1.0; //par[ parameters.Index("alphascalingC2")]; // MATCH THIS IN IMPACTFACTOR ALPHA_S WHEN NLO
    anomalous_dimension = 1.0; //par[ parameters.Index("anomalous_dimension")];
    icx0_bk = 0.01; //par[ parameters.Index("initialconditionX0")];
    sigma02 = 1.0;
    }
    double e_c          = ic_ec; //par[ parameters.Index("e_c")];
    double icx0_nlo_impfac = icx0_if; //par[ parameters.Index("initialconditionX0")];
    double initialconditionY0  = 0; //par[ parameters.Index("initialconditionY0")];
    double icTypicalPartonVirtualityQ0sqr  = 1.0; //par[ parameters.Index("icTypicalPartonVirtualityQ0sqr")];
    bool useMasses = nlodis_config::USE_MASSES;
    double qMass_light  = 0.14; // GeV --- doesn't improve fit at LO
    double qMass_u = 0.00216; // GeV, literature value
    double qMass_d = 0.00467; // GeV, literature value
    double qMass_s = 0.093; // GeV, literature value
    double qMass_charm = 1.27;
    // double qMass_b = 4.180; // GeV, literature value
    double qMass_b = 4.75; // GeV, suggested pole mass scheme standard value from Heikki

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
         << ", q_mass=" << qmass
         << endl;
    /*
    // ***Solve BK***
    */

    // double maxy = 15.4249; // from the paper, fcBK_MV.dat, 15=10+5 extra for z2improved extended evolution

    // double maxy = std::log(icx0_bk/(1e-5)) + initialconditionY0; // divisor=smallest HERA xbj in Q^2 range (1E-05)?
    // if (useImprovedZ2Bound){maxy += 5;}

    // double maxy = std::log(icx0_bk/(1e-7)); // from the paper, fcBK_MV.dat, 15=10+5 extra for z2improved extended evolution
    // if (useImprovedZ2Bound){maxy += std::log(100/1);} // extra evolution from the term log(Q^2 / Q_0^2)

    double eta0 = 0;

    double min_xbj = 1e-5;
    // NLODISSIGMAR IC + SOLVER CODE
    MV ic;                                            // Initial condition
    ic.SetQsqr(qs0sqr);
    ic.SetAnomalousDimension(anomalous_dimension);
    ic.SetE(e_c);                                     // e_c of MVe parametrization

    Dipole dipole(&ic);
    dipole.SetX0(icx0_bk);
    BKSolver solver(&dipole);
    double maxy = std::log(icx0_bk/min_xbj) + initialconditionY0; // divisor=smallest HERA xbj in Q^2 range (1E-05)?
    
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
    // if (use_custom_prescription == true){
    //     // DEFINE CUSTOM PRESCRIPTION HERE BY OVERWRITING SOME EFFECTS OF MetaPrescriptionSetter
    //     cout << "############## USE_CUSTOM_PRESCRIPTION IS SET TO TRUE ################" << endl;

    //     if (custom_presc == "cust1"){
    //         /* z2sim vs z2imp + kinematical rapidity shift comparison
    //         *  These settings force rapidity shift on with any BK evolution
    //         */
    //         cout << "# custom1 -- z2sim vs z2imp+rho: forcing eta rapidity shift (rho) to be on. USE WITH Z2IMP." << endl;

    //         ComputeSigmaR::xrapidity_funpointer x_lo_y_eta_rap_ptr; // only sub scheme LO term should have kinematical rapidity shift. No shift with unsub!
    //         ComputeSigmaR::xrapidity_NLO_funpointer x_nlo_fun_ptr;
    //         if (nlodis_config::SUB_SCHEME == nlodis_config::SUBTRACTED){
    //             // sub scheme has evolution in the LO term so it needs the kinematical shift there as well.
    //             x_lo_y_eta_rap_ptr = &ComputeSigmaR::Xrpdty_LO_targetETA;
    //             cout << "# Using shifted target ETA rapidity in SUB leading order term." << endl;
    //         }else{
    //             // unsub scheme has no evolution in the lowest order term, so no shift there, use the Y rapidity technology that doesn't shift.
    //             // Initial condition will be defined in x_eta = x0_bk.
    //             x_lo_y_eta_rap_ptr = &ComputeSigmaR::Xrpdty_LO_projectileY;
    //             cout << "# Using unshifted target ETA rapidity in UNSUB lowest order term." << endl;
    //         }
    //         x_nlo_fun_ptr = &ComputeSigmaR::Xrpdty_NLO_targetETA;
    //         cout << "# Using target ETA evolution rapidity in NLO terms" << endl;
        
    //         SigmaComputer.SetEvolutionX_LO(x_lo_y_eta_rap_ptr);
    //         SigmaComputer.SetEvolutionX_DIP(x_lo_y_eta_rap_ptr); // dipole term is always evaluated at the same rapidity as the lowest order term.
    //         SigmaComputer.SetEvolutionX_NLO(x_nlo_fun_ptr);

    //         // trbk rho prescription
    //         SigmaComputer.SetTRBKRhoPrescription(nlodis_config::TRBK_RHO_QQ0);
    //         cout << "# Using rapidity shift RHO with ANY evolution equation." << endl;
    //     }
    //     else{
    //         cout << "Unknown custom prescription: " << custom_presc << endl;
    //     }

    // }
    
    // CUBA Monte Carlo integration library algorithm setter
    SigmaComputer.SetCubaMethod(cubaMethod);

    cout << "# === Computing Reduced Cross sections ===" << endl;

            // Print column titles.
            if(true){
            #pragma omp critical
            cout    << setw(15) << "# xbj"        << " "
                    << setw(15) << "Q^2"          << " "
                    << setw(15) << "FL_data"       << " "
                    << setw(15) << "FL_err"   << " "
                    << setw(15) << "FL_theo"        << " "
                    // << setw(15) << "FL_IC"        << " "
                    // << setw(15) << "FL_LO"        << " "
                    // << setw(15) << "FL_dip"       << " "
                    // << setw(15) << "FL_qg"        << " "
                    // << setw(15) << "FL_sigma3"    << " "
                    // << setw(15) << "FT_IC"        << " "
                    // << setw(15) << "FT_LO"        << " "
                    // << setw(15) << "FT_dip"       << " "
                    // << setw(15) << "FT_qg"        << " "
                    // << setw(15) << "FT_sigma3"    << " "
                    << endl;
                    }



    // #ifdef PARALLEL_CHISQR
    //     #pragma omp parallel for schedule(dynamic) reduction(+:Q) reduction(+:xbj)
    // #endif

            // for (int i=0; i<datasets[dataset]->NumOfPoints(); i++)


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
            // #pragma omp critical
            // cout << "\r" << i+1 << "/" << totalpoints << flush;

            // Index for this in the final data array
            int dataind=0;
            for (int dseti=0; dseti < dataset; dseti++)
                dataind += datasets[dseti]->NumOfPoints();
            dataind += i;
            
            
            
            double xbj      = datasets[dataset]->xbj(i);
            double Q2       = datasets[dataset]->Qsqr(i);
            double Q        = sqrt(Q2);
            double FL_data   = datasets[dataset]->ReducedCrossSection(i);
            double FL_err = datasets[dataset]->ReducedCrossSectionError(i);
            
            datavals[dataind] = FL_data;
            dataerrs[dataind] = FL_err;
            var_xbj[dataind] = xbj;
            var_qsqr[dataind] = Q2;

            double FL=0;
            double FL_IC=0, FL_LO=0, FL_dip=0, FL_qg=0, FL_sigma3=0;
            double FT_IC=0, FT_LO=0, FT_dip=0, FT_qg=0, FT_sigma3=0;
            int calccount=0;
            if (!computeNLO && !useMasses) // Compute reduced cross section using leading order impact factors
            {
                FL_LO = SigmaComputer.Structf_LLO(Q,xbj);
                // FT_LO = SigmaComputer.Structf_TLO(Q,xbj);
                FL = FL_LO;
                ++calccount;
            }
            if (!computeNLO && useMasses)
            {
                FL_LO = SigmaComputer.Structf_LLO_massive(Q, xbj, qmass);
                FL = FL_LO;
                ++calccount;
            }

            if (computeNLO && !useSUB) // UNSUB SCHEME Full NLO impact factors for reduced cross section
            {
                if (useMasses){
                    if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
                        if (nlodis_config::MASS_SCHEME == nlodis_config::LIGHT_PLUS_CHARM_AND_BEAUTY){
                            FL_IC = SigmaComputer.Structf_LLO(Q,icx0_bk) + SigmaComputer.Structf_LLO_massive(Q,icx0_bk, qmass) + SigmaComputer.Structf_LLO_massive(Q,icx0_bk, qMass_b);
                            // FL_LO = SigmaComputer.Structf_LLO_massive(Q, xbj, qmass);
                            FL_dip = SigmaComputer.Structf_LNLOdip(Q, xbj) + SigmaComputer.Structf_LNLOdip_massive(Q, xbj, qmass) + SigmaComputer.Structf_LNLOdip_massive(Q, xbj, qMass_b);
                            FL_qg  = SigmaComputer.Structf_LNLOqg_unsub(Q, xbj) + SigmaComputer.Structf_LNLOqg_unsub_massive(Q, xbj, qmass) + SigmaComputer.Structf_LNLOqg_unsub_massive(Q, xbj, qMass_b);
                            FL = FL_IC + FL_dip + FL_qg;
                            ++calccount;
                        }
                    }else{
                        cout << "Bound loop not supported for heavy quarks!. Exit." << endl;
                        exit(1);
                    }
                }
                if (!useMasses){
                    if (useBoundLoop){
                        FL_IC = SigmaComputer.Structf_LLO(Q,icx0_bk);
                        // FT_IC = SigmaComputer.Structf_TLO(Q,icx0_bk);
                        FL_LO = SigmaComputer.Structf_LLO(Q,xbj);
                        // FT_LO = SigmaComputer.Structf_TLO(Q,xbj);
                        FL_dip = SigmaComputer.Structf_LNLOdip_z2(Q,xbj);
                        // FT_dip = SigmaComputer.Structf_TNLOdip_z2(Q,xbj);
                        FL_qg  = SigmaComputer.Structf_LNLOqg_unsub(Q,xbj);
                        // FT_qg  = SigmaComputer.Structf_TNLOqg_unsub(Q,xbj);
                        FL = FL_IC + FL_dip + FL_qg;
                        ++calccount;}
                    if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
                        FL_IC = SigmaComputer.Structf_LLO(Q,icx0_bk);
                        // FT_IC = SigmaComputer.Structf_TLO(Q,icx0_bk);
                        FL_LO = SigmaComputer.Structf_LLO(Q,xbj);
                        // FT_LO = SigmaComputer.Structf_TLO(Q,xbj);
                        FL_dip = SigmaComputer.Structf_LNLOdip(Q,xbj);
                        // FT_dip = SigmaComputer.Structf_TNLOdip(Q,xbj);
                        FL_qg  = SigmaComputer.Structf_LNLOqg_unsub(Q,xbj);
                        // FT_qg  = SigmaComputer.Structf_TNLOqg_unsub(Q,xbj);
                        FL = FL_IC + FL_dip + FL_qg;
                        ++calccount;}
                    if (useSigma3){
                        FL_sigma3 = SigmaComputer.Structf_LNLOsigma3(Q,xbj);
                        // FT_sigma3 = SigmaComputer.Structf_TNLOsigma3(Q,xbj);
                        }
                }
            }

            if (computeNLO && useSUB) // SUB SCHEME Full NLO impact factors for reduced cross section
            {
                if (useBoundLoop){
                    FL_LO = SigmaComputer.Structf_LLO(Q,xbj);
                    // FT_LO = SigmaComputer.Structf_TLO(Q,xbj);
                    FL_dip = SigmaComputer.Structf_LNLOdip_z2(Q,xbj);
                    // FT_dip = SigmaComputer.Structf_TNLOdip_z2(Q,xbj);
                    FL_qg  = SigmaComputer.Structf_LNLOqg_sub(Q,xbj);
                    // FT_qg  = SigmaComputer.Structf_TNLOqg_sub(Q,xbj);
                    FL = FL_LO + FL_dip + FL_qg;
                    ++calccount;}
                if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
                    FL_LO = SigmaComputer.Structf_LLO(Q,xbj);
                    // FT_LO = SigmaComputer.Structf_TLO(Q,xbj);
                    FL_dip = SigmaComputer.Structf_LNLOdip(Q,xbj);
                    // FT_dip = SigmaComputer.Structf_TNLOdip(Q,xbj);
                    FL_qg  = SigmaComputer.Structf_LNLOqg_sub(Q,xbj);
                    // FT_qg  = SigmaComputer.Structf_TNLOqg_sub(Q,xbj);
                    FL = FL_LO + FL_dip + FL_qg;
                    ++calccount;}
            }

            if (calccount>1)
            {
                cerr << "ERROR: Multiple computations. abort." << "count="<< calccount << endl;
                exit(1);
            }

            // Output for plotting
            if(true){
            #pragma omp critical
            cout    << setw(15) << xbj                  << " "
                    << setw(15) << Q*Q                  << " "
                    << setw(15) << FL_data              << " "
                    << setw(15) << FL_err               << " "
                    << setw(15) << sigma02*FL           << " "
                    // << setw(15) << sigma02*FL_IC        << " "
                    // << setw(15) << sigma02*FL_LO        << " "
                    // << setw(15) << sigma02*FL_dip       << " "
                    // << setw(15) << sigma02*FL_qg        << " "
                    // << setw(15) << sigma02*FL_sigma3    << " "
                    // << setw(15) << sigma02*FT_IC        << " "
                    // << setw(15) << sigma02*FT_LO        << " "
                    // << setw(15) << sigma02*FT_dip       << " "
                    // << setw(15) << sigma02*FT_qg        << " "
                    // << setw(15) << sigma02*FT_sigma3    << " "
                    << endl;
                    }
        }
        // cout << "End was reached." << endl;
        return 0;
    }
}