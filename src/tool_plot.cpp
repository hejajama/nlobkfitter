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
    // if (gsl_errno == 15 or gsl_errno == 16) return;
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
        nlodis_config::CUBA_EPSREL = 10e-3;
        nlodis_config::CUBA_MAXEVAL= 1e7;
        nlodis_config::MINR = 1e-5;
        nlodis_config::MAXR = 50;
        nlodis_config::PRINTDATA = true;
        bool useNLO = true;
        bool computeNLO = useNLO;
        string cubaMethod = "suave";

        // Oliko paperissa LO + fc BK? Taisi olla ja resummaukset vasta olivat kiinnostuksen alla sen jälkeen?
        config::LO_BK = true;  // Solve LO BK with running coupling, overrides RESUM settings
        
        config::VERBOSE = true;
        config::RINTPOINTS = 512;
        config::THETAINTPOINTS = 512;

        config::INTACCURACY = 100e-3;//0.02;
        config::MCINTACCURACY = 100e-3;//0.02;
        config::MCINTPOINTS = 1e7;
        config::MINR = 1e-5;
        config::MAXR = 50;
        config::RPOINTS = 100;
        config::DE_SOLVER_STEP = 0.4; // Rungekutta step

        // Constants
        config::NF=3;   // Only light quarks
        config::LAMBDAQCD = 0.241;


    // Data data;
    // data.SetMinQsqr(0.75);
    // data.SetMaxQsqr(50);
    // data.SetMaxX(0.01);
    // data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);

    MnUserParameters parameters;

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
            config::RESUM_DLOG = true; // Resum doulbe logs
            config::RESUM_SINGLE_LOG = true; // Resum single logs
            config::KSUB = 0.65;  // Optimal value for K_sub
            config::NO_K2 = true;  // Do not include numerically demanding full NLO part
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

/*
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
*/

/*
This needs to go as well.

    NLODISFitter fitter(parameters);
    fitter.AddDataset(data);
    fitter.SetNLO(useNLO);
    fitter.SetSUB(useSUB);
    fitter.SetSigma3(useSigma3);
    fitter.UseImprovedZ2Bound(useImprovedZ2Bound);
    fitter.UseConsistentlyBoundLoopTerm(useBoundLoop);
    fitter.SetCubaMethod(cubaMethod);
*/

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
    // cout << "=== Initial parameters ===" << endl;
    // cout << parameters << endl;
    // if(nlodis_config::VERBOSE) cout << "=== Starting fit ===" << endl;

    // params to IC
    double qs0sqr       = 0.2; //par[ parameters.Index("qs0sqr")];
    double e_c          = 1.0; //par[ parameters.Index("e_c")];
    //double fitsigma0    = 2.568*par[ parameters.Index("fitsigma0")];  // 1mb = 2.568 GeV² -- unit change into GeV
    double alphas_scaling     = 1.0; //par[ parameters.Index("alphascalingC2")]; // MATCH THIS IN IMPACTFACTOR ALPHA_S WHEN NLO
    double anomalous_dimension = 1.0; //par[ parameters.Index("anomalous_dimension")];
    double initialconditionX0  = 0.01; //par[ parameters.Index("initialconditionX0")];
    double initialconditionY0  = 0; //par[ parameters.Index("initialconditionY0")];
    double icTypicalPartonVirtualityQ0sqr  = 1.0; //par[ parameters.Index("icTypicalPartonVirtualityQ0sqr")];
    double qMass_light  = 0.14; // GeV --- doesn't improve fit at LO
    double qMass_charm = 1.35;
    bool useMasses = true;
    bool useCharm = false;

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
    // double maxy = std::log(initialconditionX0/(1e-5)) + initialconditionY0; // divisor=smallest HERA xbj in Q^2 range (1E-05)?
    double maxy = 10; // from the paper, fcBK_MV.dat

    double eta0 = 0;
    solver.SetAlphasScaling(alphas_scaling);
    solver.SetEta0(eta0);
    solver.Solve(maxy);                                // Solve up to maxy

    // Give solution to the AmplitudeLib object
    AmplitudeLib DipoleAmplitude(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());
    DipoleAmplitude.SetX0(initialconditionX0);
    DipoleAmplitude.SetOutOfRangeErrors(false);
    AmplitudeLib *DipolePointer = &DipoleAmplitude;

    ComputeSigmaR SigmaComputer(DipolePointer);
    SigmaComputer.SetX0(initialconditionX0);
    SigmaComputer.SetY0(initialconditionY0);
    SigmaComputer.SetQ0Sqr(icTypicalPartonVirtualityQ0sqr);
    SigmaComputer.SetQuarkMassLight(qMass_light);
    SigmaComputer.SetQuarkMassCharm(qMass_charm);
    // NLO: set runningcoupling and C2=Csq for the object.
    ComputeSigmaR::CmptrMemFn alphas_temppointer;
    ComputeSigmaR::CmptrMemFn_void alphas_temppointer_QG;
    if      (nlodis_config::RC_DIS == nlodis_config::DIS_RC_FIXED){
        alphas_temppointer = &ComputeSigmaR::alpha_bar_fixed;
        alphas_temppointer_QG  = &ComputeSigmaR::alpha_bar_QG_fixed;
        cout << "Using FIXED_LO" << endl;}
    else if (nlodis_config::RC_DIS == nlodis_config::DIS_RC_PARENT){
        alphas_temppointer = &ComputeSigmaR::alpha_bar_running_pd;
        alphas_temppointer_QG  = &ComputeSigmaR::alpha_bar_QG_running_pd;
        cout << "Using parent dipole RC" << endl;}
    else if (nlodis_config::RC_DIS == nlodis_config::DIS_RC_GUILLAUME){
        alphas_temppointer = &ComputeSigmaR::alpha_bar_running_pd;
        alphas_temppointer_QG  = &ComputeSigmaR::alpha_bar_QG_running_guillaume;
        cout << "Using Guillaume RC" << endl;}
    else {cout << "Problem with the choice of runnincoupling. Unkonwn config::RC_LO." << endl;}
    SigmaComputer.SetRunningCoupling(alphas_temppointer);
    SigmaComputer.SetRunningCoupling_QG(alphas_temppointer_QG);
    SigmaComputer.SetAlphasScalingC2(alphas_scaling);
    // Set z2 lower limit settings
    ComputeSigmaR::z2funpointer z2bound_funptr;
    if  (useImprovedZ2Bound == true) {z2bound_funptr = &ComputeSigmaR::z2bound_improved;}
    else                             {z2bound_funptr = &ComputeSigmaR::z2bound_simple;}
    SigmaComputer.SetImprovedZ2Bound(z2bound_funptr);
    SigmaComputer.SetSigma3BKKernel(&ComputeSigmaR::K_resum);
    SigmaComputer.SetCubaMethod(cubaMethod);


    cout << "=== Computing Reduced Cross sections ===" << endl;

            // Print column titles.
            if(true){
            #pragma omp critical
            cout    << setw(15) << "xbj"          << " "
                    << setw(15) << "Q^2"          << " "
                    << setw(15) << "FL_LO"        << " "
                    << setw(15) << "FL_dip"       << " "
                    << setw(15) << "FL_qg"        << " "
                    << setw(15) << "FL_sigma3"    << " "
                    << setw(15) << "FT_LO"        << " "
                    << setw(15) << "FT_dip"       << " "
                    << setw(15) << "FT_qg"        << " "
                    << setw(15) << "FT_sigma3"    << " "
                    << endl;
                  }


    /*
     * Loop over datapoints and compute theoretical predictions
     */
    // int points=0, totalpoints = 0;
    // for (unsigned int dataset=0; dataset<datasets.size(); dataset++)
    //     totalpoints += datasets[dataset]->NumOfPoints();

    // These loops are trivially parallerizable
    // We only parallerize the inner loop where we have about
    // 250 points (total sigmar) and 50 points (charm)
    // for (unsigned int dataset=0; dataset<datasets.size(); dataset++)
    // {
// #ifdef PARALLEL_CHISQR
//     #pragma omp parallel for schedule(dynamic) reduction(+:Q) reduction(+:xbj)
// #endif

        // for (int i=0; i<datasets[dataset]->NumOfPoints(); i++)


// Loop over grid here.
/*
 *  In the paper we used the following grids for the plots   
 *  	- Q[0] = 1.0, Q*=10^(1/20), while Q <= 10
 *		- xbj[0] = 1e-3, xbj*=10^(1/4), while xbj >= 5e-7
 * 
 */
    double Q = 1.0;
    double icx0 = 1e-2;
    double xbj = icx0;

    // #pragma omp parallel for collapse(2)
    // for (double Q = 1.0; Q <= 10 + 1e-3; Q*=pow(10,1/20) )
    for (int i=0; i<=20; i++)
    {
        // for (double xbj = 1e-3; xbj >= 5e-7 - 1e-9; xbj*=pow(10,1/4) )
        for (int j=0; j<=17; j++)
        {
            Q = 1.0*pow(10,(double)i/20.0);
            xbj = icx0/pow(10,(double)j/4.0);
            #pragma omp critical
            cout << "Q=" << Q << ", xbj=" << xbj << endl;
            
            double FL_LO=0, FL_dip=0, FL_qg=0, FL_sigma3=0;
            double FT_LO=0, FT_dip=0, FT_qg=0, FT_sigma3=0;
            int calccount=0;
            if (!computeNLO && !useMasses) // Compute reduced cross section using leading order impact factors
            {
                // theory = (fitsigma0)*SigmaComputer.SigmarLO(Q , xbj , y );
                FL_LO = SigmaComputer.Structf_LLO(Q,xbj);
                FT_LO = SigmaComputer.Structf_TLO(Q,xbj);
                ++calccount;
            }
            if (!computeNLO && useMasses)
            {
                // theory = (fitsigma0)*SigmaComputer.SigmarLOmass(Q , xbj , y );
                double alphaem=1.0/137.0;
                double structurefunfac=1./(Sq(2*M_PI)*alphaem);
                double fac = structurefunfac*Sq(Q);
                FL_LO = fac*SigmaComputer.LLOpMass(Q,xbj,useCharm);
                FT_LO = fac*SigmaComputer.TLOpMass(Q,xbj,useCharm);
                ++calccount;
            }

            if (computeNLO && !useSUB) // UNSUB SCHEME Full NLO impact factors for reduced cross section
            {
                if (useBoundLoop){
                    // theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub_UniformZ2Bound(Q , xbj , y );
                    FL_LO = SigmaComputer.Structf_LLO(Q,initialconditionX0);
                    FT_LO = SigmaComputer.Structf_TLO(Q,initialconditionX0);
                    FL_dip = SigmaComputer.Structf_LNLOdip_z2(Q,xbj);
                    FT_dip = SigmaComputer.Structf_TNLOdip_z2(Q,xbj);
                    FL_qg  = SigmaComputer.Structf_LNLOqg_unsub(Q,xbj);
                    FT_qg  = SigmaComputer.Structf_TNLOqg_unsub(Q,xbj);
                    ++calccount;}
                if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
                    // theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub(Q , xbj , y );
                    //theory = (fitsigma0)*SigmaComputer.SigmarNLOsubRisto(Q , xbj , y );
                    FL_LO = SigmaComputer.Structf_LLO(Q,initialconditionX0);
                    FT_LO = SigmaComputer.Structf_TLO(Q,initialconditionX0);
                    FL_dip = SigmaComputer.Structf_LNLOdip(Q,xbj);
                    FT_dip = SigmaComputer.Structf_TNLOdip(Q,xbj);
                    FL_qg  = SigmaComputer.Structf_LNLOqg_unsub(Q,xbj);
                    FT_qg  = SigmaComputer.Structf_TNLOqg_unsub(Q,xbj);
                    ++calccount;}
                if (useSigma3){
                    // theory += (fitsigma0)*SigmaComputer.SigmarNLOunsub_sigma3(Q , xbj , y );
                    FL_sigma3 = SigmaComputer.Structf_LNLOsigma3(Q,xbj);
                    FT_sigma3 = SigmaComputer.Structf_TNLOsigma3(Q,xbj);
                    }
            }

            if (computeNLO && useSUB) // SUB SCHEME Full NLO impact factors for reduced cross section
            {
                if (useBoundLoop){
                    // theory = (fitsigma0)*SigmaComputer.SigmarNLOsub_UniformZ2Bound(Q , xbj , y );
                    FL_LO = SigmaComputer.Structf_LLO(Q,xbj);
                    FT_LO = SigmaComputer.Structf_TLO(Q,xbj);
                    FL_dip = SigmaComputer.Structf_LNLOdip_z2(Q,xbj);
                    FT_dip = SigmaComputer.Structf_TNLOdip_z2(Q,xbj);
                    FL_qg  = SigmaComputer.Structf_LNLOqg_sub(Q,xbj);
                    FT_qg  = SigmaComputer.Structf_TNLOqg_sub(Q,xbj);
                    ++calccount;}
                if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
                    // theory = (fitsigma0)*SigmaComputer.SigmarNLOsub(Q , xbj , y );
                    //theory = (fitsigma0)*SigmaComputer.SigmarNLOsubRisto(Q , xbj , y );
                    FL_LO = SigmaComputer.Structf_LLO(Q,xbj);
                    FT_LO = SigmaComputer.Structf_TLO(Q,xbj);
                    FL_dip = SigmaComputer.Structf_LNLOdip(Q,xbj);
                    FT_dip = SigmaComputer.Structf_TNLOdip(Q,xbj);
                    FL_qg  = SigmaComputer.Structf_LNLOqg_sub(Q,xbj);
                    FT_qg  = SigmaComputer.Structf_TNLOqg_sub(Q,xbj);

                    ++calccount;}
            }

            if (calccount>1)
            {
              cerr << "ERROR: Multiple computations. abort." << "count="<< calccount << endl;
              exit(1);
            }

            // if (std::isnan(theory) or std::isinf(theory))
            // {
            //     cerr << "Warning: theory result " << theory << " with parameters " << PrintVector(par) << endl;
            //     exit(1);
            // }

            // chisqr += datasets[dataset]->Weight()*SQR( (theory - sigmar) / sigmar_err );
            // points = points + datasets[dataset]->Weight();

            // Output for plotting
            if(true){
            #pragma omp critical
            cout    << setw(15) << xbj          << " "
                    << setw(15) << Q*Q          << " "
                    << setw(15) << FL_LO        << " "
                    << setw(15) << FL_dip       << " "
                    << setw(15) << FL_qg        << " "
                    << setw(15) << FL_sigma3    << " "
                    << setw(15) << FT_LO        << " "
                    << setw(15) << FT_dip       << " "
                    << setw(15) << FT_qg        << " "
                    << setw(15) << FT_sigma3    << " "
                    << endl;
                  }
        }
    }
    // }
    // cout << endl << "# Calculated chi^2/N = " << chisqr/points << " (N=" << points << "), parameters (" << PrintVector(par) << ")" << endl<<endl;


    return 0;
}
