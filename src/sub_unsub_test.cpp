#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <csignal>
#include <ctime>
#include <chrono>
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
    // if (gsl_errno == 15 or gsl_errno == 16) return;
    // Ugly hack, comes from the edges of the z integral in virtual_photon.cpp
    // Overflows come from IPsat::bint when it is done analytically
    // Hope is that these errors are handled correctly everywhere
    errors_counter++;
    std::cerr << file << ":"<< line <<": Error " << errors_counter << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}

int main( int argc, char* argv[] )
// int main()
{
    gsl_set_error_handler(&ErrHandlerCustom);

    // NLO DIS SIGMA_R COMPUTATION CONFIGS
    nlodis_config::CUBA_EPSREL = 10e-3;
    nlodis_config::CUBA_MAXEVAL= 5e7;
    nlodis_config::MINR = 1e-6;
    nlodis_config::MAXR = 40;
    nlodis_config::PRINTDATA = true;
    bool useNLO = true;
    bool computeNLO = useNLO;
    // string cubaMethod = "vegas";
    string cubaMethod = "suave";

    config::NO_K2 = true;  // Do not include numerically demanding full NLO part
    config::KINEMATICAL_CONSTRAINT = config::KC_NONE;

    config::VERBOSE = true;
    //config::RINTPOINTS = 512/4;
    //config::THETAINTPOINTS = 512/4;

    config::INTACCURACY = 10e-3;//0.02;
    config::MINR = 1e-6;
    config::MAXR = 40;
    config::RPOINTS = 100;
    config::DE_SOLVER_STEP = 0.4; // Rungekutta step

    // Constants
    config::NF=3;   // Only light quarks
    config::LAMBDAQCD = 0.241;

    MnUserParameters parameters;

    bool useSUB, useResumBK, useKCBK, useImprovedZ2Bound, useBoundLoop;
    bool useSigma3 = false;
    string helpstring = "Argument order: BK RC useImprovedZ2Bound useBoundLoop Q C^2 X0 gamma Q0sq Y0 eta0\nresumbk/trbk/lobk parentrc/guillaumerc/fixedrc z2improved/z2simple z2boundloop/unboundloop";
    string string_sub, string_bk, string_rc;
    if (argc<2){ cout << helpstring << endl; return 0;}
    // Argv[0] is the name of the program

    string_bk = string(argv [1]);
    if (string_bk == "resumbk"){
            config::EULER_METHOD = false;    // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = true;       // Resum doulbe logs
            config::RESUM_SINGLE_LOG = true; // Resum single logs
            config::KSUB = 0.65;             // Optimal value for K_sub
            config::NO_K2 = true;            // Do not include numerically demanding full NLO part
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_RESUM;
    }else if (string_bk == "kcbk"){
            config::EULER_METHOD = true;     // Kinematical constraint requires this
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_BEUF_K_PLUS;
            config::DE_SOLVER_STEP = 0.05;  //0.02; // Euler method requires smaller step than RungeKutta!
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_KCBK_BEUF;
    }else if (string_bk == "trbk"){  // Target Rapidity BK
            config::EULER_METHOD = true;     // Kinematical constraint requires this
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_EDMOND_K_MINUS;
            config::DE_SOLVER_STEP = 0.05;  //0.02; // Euler method requires smaller step than RungeKutta!
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_TRBK_EDMOND;
            nlodis_config::TRBK_RHO_PRESC = nlodis_config::TRBK_RHO_MAX_X_Y_R;
    }else if (string_bk == "lobk"){
            config::EULER_METHOD = false;   // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_NONE;
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_LOBK_EXPLICIT;
    }else if (string_bk == "lobkold"){
            config::EULER_METHOD = false;   // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_NONE;
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_LOBK_Z2TOZERO;
    } else {cout << helpstring << endl; return -1;}

    string_rc = string(argv [2]);
    if (string_rc == "parentrc" or string_rc == "pdrc"){
            config::RC_LO = config::PARENT_LO;
            config::RESUM_RC = config::RESUM_RC_PARENT;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_PARENT;
    } else if (string_rc == "guillaumerc" or string_rc == "gbrc"){
            config::RC_LO = config::GUILLAUME_LO;
            config::RESUM_RC = config::RESUM_RC_GUILLAUME;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_GUILLAUME;
    } else if (string_rc == "fixedrc" or string_rc == "fc"){
            config::RC_LO = config::FIXED_LO;
            config::RESUM_RC = config::RESUM_RC_FIXED;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_FIXED;
    } else if (string_rc == "smallestrc" or string_rc == "sdrc"){
            config::RC_LO = config::SMALLEST_LO;
            config::RESUM_RC = config::RESUM_RC_BALITSKY;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_SMALLEST;
    } else {cout << helpstring << endl; return -1;}

    if (string(argv [3]) == "z2improved" or string(argv [3]) == "z2imp"){
        nlodis_config::Z2MINIMUM = nlodis_config::Z2IMPROVED;
        useImprovedZ2Bound = true;
    } else if (string(argv [3]) == "z2simple" or string(argv [3]) == "z2sim"){
        nlodis_config::Z2MINIMUM = nlodis_config::Z2SIMPLE;
        useImprovedZ2Bound = false;
    } else {cout << helpstring << endl; return -1;}

    if (string(argv [4]) == "z2boundloop" or string(argv [4]) == "z2b"){
        useBoundLoop = true;
    } else if (string(argv [4]) == "unboundloop" or string(argv [4]) == "unb"){
        useBoundLoop = false;
    } else {cout << helpstring << endl; return -1;}


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
            << "# config::INTACCURACY = " << config::INTACCURACY
                << ", config::RPOINTS = " << config::RPOINTS
                << ", config::DE_SOLVER_STEP = " << config::DE_SOLVER_STEP
                << ", config::{MINR, MAXR}, nlodis_config::{MINR, MAXR} = " << config::MINR << " " << config::MAXR << " " << nlodis_config::MINR << " " << nlodis_config::MAXR
                << endl;

/*
    *   PARAMS TO IC
    */
    double qs0sqr       = 0.2; //par[ parameters.Index("qs0sqr")];
    double alphas_scaling     = 1.0; //par[ parameters.Index("alphascalingC2")]; // MATCH THIS IN IMPACTFACTOR ALPHA_S WHEN NLO
    double anomalous_dimension = 1.0; //par[ parameters.Index("anomalous_dimension")];
    double e_c          = 1.0; //par[ parameters.Index("e_c")];
    double icx0_nlo_impfac = 1.0; //par[ parameters.Index("initialconditionX0")];
    double icx0_bk = 1.0; //par[ parameters.Index("initialconditionX0")];
    double initialconditionY0  = 0; //par[ parameters.Index("initialconditionY0")];
    double icTypicalPartonVirtualityQ0sqr  = 1.0; //par[ parameters.Index("icTypicalPartonVirtualityQ0sqr")];
    double qMass_light  = 0.14; // GeV --- doesn't improve fit at LO
    double qMass_charm = 1.35;
    bool useMasses = true;
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
         << endl;

    /*
    // ***Solve resummed BK***
    */

    // /*
    MV ic;                                            // Initial condition
    ic.SetQsqr(qs0sqr);
    ic.SetAnomalousDimension(anomalous_dimension);
    ic.SetE(e_c);                                     // e_c of MVe parametrization

    Dipole dipole(&ic);
    BKSolver solver(&dipole);
    // double maxy = std::log(initialconditionX0/(1e-5)) + initialconditionY0; // divisor=smallest HERA xbj in Q^2 range (1E-05)?
    double maxy = 10; // from the paper, fcBK_MV.dat, 15=10+5 extra for z2improved extended evolution
    if (useImprovedZ2Bound){maxy += 5;}

    double eta0 = 0;
    solver.SetAlphasScaling(alphas_scaling);
    // solver.SetEta0(eta0);
    solver.Solve(maxy);                                // Solve up to maxy

    // solver.GetDipole()->Save("./out/dipoles/dipole_lobk_fc_RK_step0.2_rpoints100_rmin1e-6rmax50_INTACC10e-3.dat");
    // cout << "Saved dipole to file, exiting." << endl;
    // exit(0);
    // */

    // Give solution to the AmplitudeLib object
    AmplitudeLib DipoleAmplitude(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());
    // AmplitudeLib DipoleAmplitude("./data/paper1dipole/pap1_fcBK_MV.dat"); // pap1_fcBK_MV.dat, pap1_rcBK_MV_parent.dat
    // AmplitudeLib DipoleAmplitude("./out/dipoles/dipole_lobk_fc_step0.2_rpoints400-2.dat"); // pap1_fcBK_MV.dat, pap1_rcBK_MV_parent.dat
    // AmplitudeLib DipoleAmplitude("./out/dipoles/dipole_lobk_fc_step0.2_rpoints400_rmin1e-6rmax50_INTACC2e-3.dat"); // pap1_fcBK_MV.dat, pap1_rcBK_MV_parent.dat
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
            if(true){
            #pragma omp critical
            cout    << setw(15) << "# xbj"               << " "
                    << setw(15) << "Q^2"             << " "
                    << setw(15) << "FL_LO_sub"         << " "
                    << setw(15) << "FL_qg_sub"         << " "
                    << setw(15) << "error_sub_L"       << " "
                    << setw(15) << "FT_LO_sub"         << " "
                    << setw(15) << "FT_qg_sub"         << " "
                    << setw(15) << "error_sub_T"       << " "
                    << setw(15) << "FL_IC_unsub"       << " "
                    << setw(15) << "FL_qg_unsub"       << " "
                    << setw(15) << "error_unsub_L"     << " "
                    << setw(15) << "FT_IC_unsub"       << " "
                    << setw(15) << "FT_qg_unsub"       << " "
                    << setw(15) << "error_unsub_T"     << " "
                    << setw(15) << "squared diff L"     << " "
                    << setw(15) << "squared diff T"     << " "
                    << endl;
                    }

    /*
    *  In the paper we used the following grids for the plots   
    *  	- Q[0] = 1.0, Q*=10^(1/20), while Q <= 10
    *		- xbj[0] = 1e-3, xbj*=10^(1/4), while xbj >= 5e-7
    * 
    */

    int numpoints=0;
    std::vector< std::tuple<int, int> > coordinates;

    // #pragma omp parallel for collapse(2)
    // for (int i=0; i<=20; i+=17)  // Q^2 = {1,50}
    for (int i=0; i<=20; i+=4)  // Q^2 in [1,100]
    // for (int i=0; i<=1; i++)
    {
        for (int j=0; j<=17; j+=4)  // xbj in [5.62341e-07, 1e-2]
        //for (int j=4; j<=12; j+=8)  // xbj = {1e-3, 1e-5}
        // for (int j=0; j<=1; j++)
        {
            if (!((i == 0 or i == 17) or (j == 4 or j == 12))) { continue; }
            coordinates.emplace_back(i,j);
            numpoints++;
        }
    }

    double icQ = 1.0;
    double icx0 = icx0_bk;
    double chisqr_L = 0;
    double chisqr_T = 0;

    float duration_L_sub;
    float duration_T_sub;
    float duration_L_unsub;
    float duration_T_unsub;
// Loop over data points and compute both SUB and UNSUB and compute some chi^2 or similar qualifier
// It suffices to compute sigma^IC, sigma^LO, sigma^qg in either scheme and compare ic+qg(unsub) and lo+qg(sub)
    #pragma omp parallel for
    for (size_t k=0; k< coordinates.size(); k++)
        {
        int i,j;
        std::tie(i,j) = coordinates[k];

        if (j==0 and cubaMethod=="suave"){continue;}
        double Q = 1.0*pow(10,(double)i/20.0);
        double xbj = icx0/pow(10,(double)j/4.0);

        double FL_LO_sub=0, FL_qg_sub=0;
        double FT_LO_sub=0, FT_qg_sub=0;
        double FL_IC_unsub=0, FL_qg_unsub=0;
        double FT_IC_unsub=0, FT_qg_unsub=0;

        double error_sub_L=0, error_unsub_L=0;
        double error_sub_T=0, error_unsub_T=0;

        // UNSUB SCHEME Full NLO impact factors for reduced cross section
        {
        auto start_l_unsub = std::chrono::high_resolution_clock::now();
        FL_IC_unsub = SigmaComputer.Structf_LLO(Q,icx0_bk);
        FL_qg_unsub  = SigmaComputer.Structf_LNLOqg_unsub(Q,xbj);
        auto stop_l_unsub = std::chrono::high_resolution_clock::now();
        auto dur_l_unsub = std::chrono::duration_cast<std::chrono::milliseconds>(stop_l_unsub - start_l_unsub);

        auto start_t_unsub = std::chrono::high_resolution_clock::now();
        FT_IC_unsub = SigmaComputer.Structf_TLO(Q,icx0_bk);
        FT_qg_unsub  = SigmaComputer.Structf_TNLOqg_unsub(Q,xbj);
        auto stop_t_unsub = std::chrono::high_resolution_clock::now();
        auto dur_t_unsub = std::chrono::duration_cast<std::chrono::milliseconds>(stop_t_unsub - start_t_unsub);

        duration_L_unsub += dur_l_unsub.count();
        duration_T_unsub += dur_t_unsub.count();

        error_unsub_L = nlodis_config::CUBA_EPSREL*sqrt( Sq(FL_IC_unsub)+Sq(FL_qg_unsub) );
        error_unsub_T = nlodis_config::CUBA_EPSREL*sqrt( Sq(FT_IC_unsub)+Sq(FT_qg_unsub) );
        }

        // SUB SCHEME Full NLO impact factors for reduced cross section
        {
        auto start_l_sub = std::chrono::high_resolution_clock::now();
        FL_LO_sub = SigmaComputer.Structf_LLO(Q,xbj);
        FL_qg_sub  = SigmaComputer.Structf_LNLOqg_sub(Q,xbj);
        auto stop_l_sub = std::chrono::high_resolution_clock::now();
        auto dur_l_sub = std::chrono::duration_cast<std::chrono::milliseconds>(stop_l_sub - start_l_sub);

        auto start_t_sub = std::chrono::high_resolution_clock::now();
        FT_LO_sub = SigmaComputer.Structf_TLO(Q,xbj);
        FT_qg_sub  = SigmaComputer.Structf_TNLOqg_sub(Q,xbj);
        auto stop_t_sub = std::chrono::high_resolution_clock::now();
        auto dur_t_sub = std::chrono::duration_cast<std::chrono::milliseconds>(stop_t_sub - start_t_sub);

        duration_L_sub += dur_l_sub.count();
        duration_T_sub += dur_t_sub.count();

        error_sub_L = nlodis_config::CUBA_EPSREL*sqrt( Sq(FL_LO_sub)+Sq(FL_qg_sub) );
        error_sub_T = nlodis_config::CUBA_EPSREL*sqrt( Sq(FT_LO_sub)+Sq(FT_qg_sub) );
        }

        // chisq = sum (theory - data)^2 / error^2
        chisqr_L += ( Sq(FL_LO_sub + FL_qg_sub - (FL_IC_unsub + FL_qg_unsub)) ) / ( Sq(error_sub_L) + Sq(error_unsub_L) );
        chisqr_T += ( Sq(FT_LO_sub + FT_qg_sub - (FT_IC_unsub + FT_qg_unsub)) ) / ( Sq(error_sub_T) + Sq(error_unsub_T) );


        //prints
        if(true){
            #pragma omp critical
            cout    << setw(15) << xbj               << " "
                    << setw(15) << Sq(Q)             << " "
                    << setw(15) << FL_LO_sub         << " "
                    << setw(15) << FL_qg_sub         << " "
                    << setw(15) << error_sub_L       << " "
                    << setw(15) << FT_LO_sub         << " "
                    << setw(15) << FT_qg_sub         << " "
                    << setw(15) << error_sub_T       << " "
                    << setw(15) << FL_IC_unsub       << " "
                    << setw(15) << FL_qg_unsub       << " "
                    << setw(15) << error_unsub_L     << " "
                    << setw(15) << FT_IC_unsub       << " "
                    << setw(15) << FT_qg_unsub       << " "
                    << setw(15) << error_unsub_T     << " "
                    << setw(15) << Sq(FL_LO_sub + FL_qg_sub - (FL_IC_unsub + FL_qg_unsub))     << " "
                    << setw(15) << Sq(FT_LO_sub + FT_qg_sub - (FT_IC_unsub + FT_qg_unsub))     << " "
                    << endl;
        }
    }

    // final prints on sub unsub agreement
    cout    << setw(15) << "# Overall SUB UNSUB agreement:"
            << endl
            << setw(15) << "# chisqr_L/N: " << chisqr_L/numpoints << " "
            << setw(15) << "# chisqr_T/N: " << chisqr_T/numpoints
            << endl
            << setw(15) << "# L sub t[s]: " << duration_L_sub/1000.0 << " "
            << setw(15) << "# T sub t[s]: " << duration_T_sub/1000.0 << " "
            << endl
            << setw(15) << "# L unsub t[s]: " << duration_L_unsub/1000.0 << " "
            << setw(15) << "# T unsub t[s]: " << duration_T_unsub/1000.0 << " "
            << endl;

}