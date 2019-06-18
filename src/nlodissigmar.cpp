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

#include "nlodis_config.hpp"
#include "nlodissigmar.hpp"
#include "helper.hpp"


namespace sigmar_config{
    // double rmax=30.0; //300.0;
    // double rmin=1e-5;
    double Nc=3.0;
    double sumef=6.0/9.0; // light quarks uds only.
    double CF=4.0/3.0; // (Nc()*Nc()-1.0)/(2.0*Nc());
    double alphaem=1.0/137.0;
    double lambdaqcd=0.241; //GeV
    double structurefunfac=1./(Sq(2*M_PI)*alphaem);
}
using namespace sigmar_config;

///===========================================================================================
//// MAIN FUNCTIONS HERE

// at low Q² Z exchange negligible ( we are only modeling gamma exchange) and the reduced cross section simplifies into
// sigma_r = F_2 - f(y) F_L = F_2 - y^2 / Y_+ * F_L = F_2 - y^2/(1+(1-y)^2) * F_L
double ComputeSigmaR::SigmarLO (double Q , double xbj, double y) {
    double FL = Structf_LLO(Q,xbj);
    double FT = Structf_TLO(Q,xbj);
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigmar = F2 - fy*FL;
    return sigmar;
}
double ComputeSigmaR::SigmarLOmass (double Q , double xbj, double y, bool charm) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LLOpMass(Q,xbj,charm);
    double FT = fac*TLOpMass(Q,xbj,charm);
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

// UNSUB ---------------------
double ComputeSigmaR::SigmarNLOunsub ( double Q , double xbj, double y) {
    double FLic  = Structf_LLO(Q,icX0);
    double FTic  = Structf_TLO(Q,icX0);
    double FLdip = Structf_LNLOdip(Q,xbj);
    double FTdip = Structf_TNLOdip(Q,xbj);
    double FLqg  = Structf_LNLOqg_unsub(Q,xbj);
    double FTqg  = Structf_TNLOqg_unsub(Q,xbj);
    double FL = FLic + FLdip + FLqg;
    double FT = FTic + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

double ComputeSigmaR::SigmarNLOunsub_UniformZ2Bound ( double Q , double xbj, double y) {
    double FLic  = Structf_LLO(Q,icX0);
    double FTic  = Structf_TLO(Q,icX0);
    double FLdip = Structf_LNLOdip_z2(Q,xbj);
    double FTdip = Structf_TNLOdip_z2(Q,xbj);
    double FLqg  = Structf_LNLOqg_unsub(Q,xbj);
    double FTqg  = Structf_TNLOqg_unsub(Q,xbj);
    double FL = FLic + FLdip + FLqg;
    double FT = FTic + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

double ComputeSigmaR::SigmarNLOunsub_sigma3 ( double Q , double xbj, double y) {
    double FL3  = Structf_LNLOsigma3(Q,xbj);
    double FT3  = Structf_TNLOsigma3(Q,xbj);
    double F2 = FL3+FT3;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL3;
    return sigma;
}


// SUB ---------------------
double ComputeSigmaR::SigmarNLOsub ( double Q , double xbj, double y) {
    double FLlo  = Structf_LLO(Q,xbj);
    double FTlo  = Structf_TLO(Q,xbj);
    double FLdip = Structf_LNLOdip(Q,xbj);
    double FTdip = Structf_TNLOdip(Q,xbj);
    double FLqg  = Structf_LNLOqg_sub(Q,xbj);
    double FTqg  = Structf_TNLOqg_sub(Q,xbj);
    double FL = FLlo + FLdip + FLqg;
    double FT = FTlo + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

double ComputeSigmaR::SigmarNLOsub_UniformZ2Bound ( double Q , double xbj, double y) {
    double FLlo  = Structf_LLO(Q,xbj);
    double FTlo  = Structf_TLO(Q,xbj);
    double FLdip = Structf_LNLOdip_z2(Q,xbj);
    double FTdip = Structf_TNLOdip_z2(Q,xbj);
    double FLqg  = Structf_LNLOqg_sub(Q,xbj);
    double FTqg  = Structf_TNLOqg_sub(Q,xbj);
    double FL = FLlo + FLdip + FLqg;
    double FT = FTlo + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

// RISTO formulations ---------------------
double ComputeSigmaR::SigmarNLOunsubRisto ( double Q , double xbj, double y) {
    double fac = structurefunfac*Sq(Q);
    double FLic  = fac*LLOp(Q,icX0);
    double FTic  = fac*TLOp(Q,icX0);
    double FLdip = fac*LNLOdip(Q,xbj);
    double FTdip = fac*TNLOdip(Q,xbj);
    double FLqg = fac*LNLOqgunsubRisto(Q,xbj);
    double FTqg = fac*TNLOqgunsubRisto(Q,xbj);
    double FL = FLic + FLdip + FLqg;
    double FT = FTic + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

double ComputeSigmaR::SigmarNLOsubRisto ( double Q , double xbj, double y) {
    double fac = structurefunfac*Sq(Q);
    double FLlo  = fac*LLOp(Q,xbj);
    double FTlo  = fac*TLOp(Q,xbj);
    double FLdip = fac*LNLOdip(Q,xbj);
    double FTdip = fac*TNLOdip(Q,xbj);
    double FLqg = fac*LNLOqgsubRisto(Q,xbj);
    double FTqg = fac*TNLOqgsubRisto(Q,xbj);
    double FL = FLlo + FLdip + FLqg;
    double FT = FTlo + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

// STRUCTURE FUNCTIONS -------------------------
double ComputeSigmaR::Structf_LLO ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LLOp( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_TLO ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TLOp( Q , xbj );
    return FT;
}

double ComputeSigmaR::Structf_LNLOdip ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LNLOdip( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_TNLOdip ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TNLOdip( Q , xbj );
    return FT;
}

double ComputeSigmaR::Structf_LNLOdip_z2 ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LNLOdip_z2( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_TNLOdip_z2 ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TNLOdip_z2( Q , xbj );
    return FT;
}

double ComputeSigmaR::Structf_LNLOqg_sub ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LNLOqgsub( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_LNLOqg_unsub ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LNLOqgunsub( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_TNLOqg_sub ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TNLOqgsub( Q , xbj );
    return FT;
}

double ComputeSigmaR::Structf_TNLOqg_unsub ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TNLOqgunsub( Q , xbj );
    return FT;
}

double ComputeSigmaR::Structf_LNLOsigma3 ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LNLOsigma3( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_TNLOsigma3 ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TNLOsigma3( Q , xbj );
    return FT;
}


// obsolete?
double ComputeSigmaR::Structf_LFULLNLOunsub ( double Q , double xbj) {
    double fac = structurefunfac*Sq(Q);
    double FLic  = fac*LLOp(Q,icX0);
    //double FLdip = fac*LNLOdip(Q,xbj);
    double FLqg = fac*LNLOqgunsub(Q,xbj);
    double FL = FLic + FLqg;
    return FL;
}

double ComputeSigmaR::Structf_TFULLNLOunsub ( double Q , double xbj) {
    double fac = structurefunfac*Sq(Q);
    double FTic  = fac*TLOp(Q,icX0);
    //double FTdip = fac*TNLOdip(Q,xbj);
    double FTqg = fac*TNLOqgunsub(Q,xbj);
    double FT = FTic + FTqg;
    return FT;
}

double ComputeSigmaR::Structf_LFULLNLOsub ( double Q , double xbj) {
    double fac = structurefunfac*Sq(Q);
    double FLlo  = fac*LLOp(Q,xbj);
    //double FLdip = fac*LNLOdip(Q,xbj);
    double FLqg = fac*LNLOqgsub(Q,xbj);
    double FL = FLlo + FLqg;
    return FL;
}

double ComputeSigmaR::Structf_TFULLNLOsub ( double Q , double xbj) {
    double fac = structurefunfac*Sq(Q);
    double FTlo  = fac*TLOp(Q,xbj);
    //double FTdip = fac*TNLOdip(Q,xbj);
    double FTqg = fac*TNLOqgsub(Q,xbj);
    double FT = FTlo + FTqg;
    return FT;
}


///===========================================================================================
double NLODISFitter::operator()(const std::vector<double>& par) const
{
    double chisqr = 0;
    // params to IC
    double qs0sqr       = par[ parameters.Index("qs0sqr")];
    double e_c          = par[ parameters.Index("e_c")];
    //double fitsigma0    = 2.568*par[ parameters.Index("fitsigma0")];  // 1mb = 2.568 GeV² -- unit change into GeV
    double alphas_scaling       = par[ parameters.Index("alphascalingC2")]; // MATCH THIS IN IMPACTFACTOR ALPHA_S WHEN NLO
    double anomalous_dimension  = par[ parameters.Index("anomalous_dimension")];
    double icx0_nlo_impfac  = par[ parameters.Index("icx0_nlo_impfac")];
    double icx0_bk  = par[ parameters.Index("icx0_bk")];
    double initialconditionY0  = par[ parameters.Index("initialconditionY0")];
    double icTypicalPartonVirtualityQ0sqr  = par[ parameters.Index("icTypicalPartonVirtualityQ0sqr")];
    double qMass_light  = 0.14; // GeV --- doesn't improve fit at LO
    double qMass_charm = 1.35;
    bool useMasses = true;
    bool useCharm = false;


    if (std::abs(initialconditionY0) > 1e-4)
    {
        cerr << "WARNING! initialcconditionY0" << initialconditionY0 << ", but values != 0 are not guaranteed to work! Talk to Heikki before using this!" << endl;
    }
    if (qs0sqr < 0.0001 or qs0sqr > 100 or alphas_scaling < 0.01 or alphas_scaling > 99999
        /*or fitsigma0 < 0.1 or fitsigma0 > 999 */
        or e_c < 1 or e_c > 9999)
    return 9999999;

    cout << "=== Initializing Chi^2 regression === "<< " parameters (" << PrintVector(par) << ")" << endl;
    // Manual limiting for parameter range as a fail safe for Minuit2 craziness:
    /*
    if (  (qs0sqr < 0)    || (alphas_scaling < 0)     || (anomalous_dimension < 0) || (e_c < 0.4) || (initialconditionX0 < 0.01) ||
            (qs0sqr > 100)  || (alphas_scaling > 1000)  || (anomalous_dimension > 2) ) {
            chisqr = 1e7 * ( 1.0 + abs(qs0sqr) + abs(alphas_scaling) + abs(anomalous_dimension) + abs(e_c - 0.4)); // Should be some crazy large chisqr to deter Minuit2 from using this kind of parametrization.
            cout << endl << "# Calculated ILL PARAMETER chi^2 = " << chisqr  << ", parameters (" << PrintVector(par) << ")" << endl<<endl;
            return chisqr;
        }
    */

    /*
    // ***Solve resummed BK***
    */
    // cout << "=== Initialize BK solver ===" << endl;
   
    MV ic;                                            // Initial condition
    ic.SetQsqr(qs0sqr);
    ic.SetAnomalousDimension(anomalous_dimension);
    ic.SetE(e_c);                                     // e_c of MVe parametrization

    Dipole dipole(&ic);
    dipole.SetX0(icx0_bk);
    BKSolver solver(&dipole);
    double maxy = std::log(icx0_bk/(1e-5)) + initialconditionY0; // divisor=smallest HERA xbj in Q^2 range (1E-05)?

    // cout << "=== Solving BK ===" << endl;

    solver.SetAlphasScaling(alphas_scaling);
    solver.SetEta0(par[ parameters.Index("eta0")]);
    solver.SetX0(icx0_bk);
    solver.SetICX0_nlo_impfac(icx0_nlo_impfac);
    solver.SetICTypicalPartonVirtualityQ0sqr(icTypicalPartonVirtualityQ0sqr);
    solver.SetTmpOutput("tmp_datafile_nokc.dat");
    solver.Solve(maxy);

    //solver.GetDipole()->Save("output_dipole_lobk_x0=10_euler_DESTEP0.02.dat");
    //solver.GetDipole()->Save("output_dipole_lobk_x0=10_euler_KINCOST_DESTEP0.08.dat");
    //solver.GetDipole()->Save("output_dipole_lobk_x0=10_rungekutta_DESTEP0.4.dat");
    
    // Give solution to the AmplitudeLib object
    AmplitudeLib DipoleAmplitude(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());
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
    else {
        cout << "ERROR: Problem with the choice of runnincoupling. Unkonwn nlodis_config::RC_DIS." << endl;
        exit(1);
    }
    SigmaComputer.SetRunningCoupling(alphas_temppointer);
    SigmaComputer.SetRunningCoupling_QG(alphas_temppointer_QG);
    SigmaComputer.SetAlphasScalingC2(alphas_scaling);
    // Set z2 lower limit settings
    ComputeSigmaR::z2funpointer z2bound_funptr;
    if  (useImprovedZ2Bound == true) {z2bound_funptr = &ComputeSigmaR::z2bound_improved;}
    else                             {z2bound_funptr = &ComputeSigmaR::z2bound_simple;}
    SigmaComputer.SetImprovedZ2Bound(z2bound_funptr);
    // Dipole evalution X, (y = log(x0/X))
    // if using 'sub' scheme with z2imp one must use the extended evolution variable for sigma_LO as well.
    ComputeSigmaR::xrapidity_funpointer x_fun_ptr;
    if (UseSub and useImprovedZ2Bound) {x_fun_ptr = &ComputeSigmaR::Xrpdty_LO_improved;}
    else        {x_fun_ptr = &ComputeSigmaR::Xrpdty_LO_simple;}
    SigmaComputer.SetEvolutionX_LO(x_fun_ptr);
    SigmaComputer.SetEvolutionX_DIP(x_fun_ptr);
    // sigma_3 BK correction
    SigmaComputer.SetSigma3BKKernel(&ComputeSigmaR::K_resum);
    // sub scheme subtraction term choice
    SigmaComputer.SetSubTermKernel(nlodis_config::SUB_TERM_KERNEL);
    // CUBA Monte Carlo integration library algorithm setter
    SigmaComputer.SetCubaMethod(cubaMethod);

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
                theory = (fitsigma0)*SigmaComputer.SigmarLO(Q , xbj , y );
                ++calccount;
            }
            if (!computeNLO && useMasses)
            {
                theory=0;
                if (datasets[dataset]->OnlyCharm(i)==false)
                {
                    theory = (fitsigma0)*SigmaComputer.SigmarLOmass(Q , xbj , y, false );
                }
                if (xbj*(1.0 + 4.0*1.35*1.35/(Q*Q)) < 0.01 and useCharm)
                {
                    theory_charm = (fitsigma0)*SigmaComputer.SigmarLOmass(Q , xbj*(1.0 + 4.0*1.35*1.35/(Q*Q)) , y, true ); // charm
                    if (datasets[dataset]->OnlyCharm(i) == true)
                        theory = theory_charm;
                    else
                        theory = theory + theory_charm;
                }
                  ++calccount;
            }
            
            if (computeNLO && !UseSub) // UNSUB SCHEME Full NLO impact factors for reduced cross section
            {
                if (useBoundLoop){
                    theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub_UniformZ2Bound(Q , xbj , y );
                    ++calccount;}
                if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
                    theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub(Q , xbj , y );
                    //theory = (fitsigma0)*SigmaComputer.SigmarNLOsubRisto(Q , xbj , y );
                    ++calccount;}
                if (UseSigma3){
                    theory += (fitsigma0)*SigmaComputer.SigmarNLOunsub_sigma3(Q , xbj , y );
            }
            }
            
            if (computeNLO && UseSub) // SUB SCHEME Full NLO impact factors for reduced cross section
            {
                if (useBoundLoop){
                    theory = (fitsigma0)*SigmaComputer.SigmarNLOsub_UniformZ2Bound(Q , xbj , y );
                    ++calccount;}
                if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
                    theory = (fitsigma0)*SigmaComputer.SigmarNLOsub(Q , xbj , y );
                    //theory = (fitsigma0)*SigmaComputer.SigmarNLOsubRisto(Q , xbj , y );
                    ++calccount;}
            }
            
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
    // Output for plotting
    if(nlodis_config::PRINTDATA){
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
    

    // delete DipolePointer;
    return chisqr_over_n * points;  // Return chisqr for Minuit error estimation.
}

void NLODISFitter::AddDataset(Data& d)
{
    datasets.push_back(&d);
}

NLODISFitter::NLODISFitter(MnUserParameters parameters_)
{
    parameters = parameters_;
    
    cubaMethod = "suave";  // Default choise for Cuba
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
    double border=1e-6, maxchisq=10, mindeviation=0.25;
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

///===========================================================================================
///===========================================================================================
// HELPERS
double ComputeSigmaR::Sr(double r, double x) {
    double Srx;//, Nrx;
    if(r<nlodis_config::MINR){
        Srx = 1.;
    }else if(r>nlodis_config::MAXR-1e-7){
        Srx = 0.;
    }else{
        // Srx = ClassScopeDipolePointer->S_y(r, log(1/x))
        // Note: icY0 controls how much evolution we have before we define that we have "initial condition"
        if (x > icX0_bk){
            Srx = ClassScopeDipolePointer->S(r, (icX0_bk*std::exp(-icY0)) ) ; //1-Nrx;
        } else {
            Srx = ClassScopeDipolePointer->S(r, (x*std::exp(-icY0)) ) ; //1-Nrx;
        }
    }
    return Srx;
}

double ComputeSigmaR::SrTripole(double x01, double x02, double x21, double x) {
    return Nc/(2*CF)*(Sr(x02,x)*Sr(x21,x) - 1/Sq(Nc)*Sr(x01,x));
}

double ComputeSigmaR::P(double z) {
    return 1.0-z+Sq(z)*0.5; // EI TÄSMÄÄ ARTIKKELIN P-funktion KANSSA. P=1/2*Partikkeli
}

double ComputeSigmaR::heaviside_theta(double x) {
    double result = 0;
    if (x >= 0) {result = 1.0;}
    return result;
}


///===========================================================================================
// RUNNING COUPLINGS
double ComputeSigmaR::alpha_bar_running_pd_sharpcutoff( double rsq ) { // alphabar = Nc/M_PI * alphas
    double logvar = gsl_sf_log(	(4*(alpha_scaling_C2_)/(rsq*Sq(lambdaqcd))	)	);
    double AlphaSres = 12*M_PI/(	(33.0-2*Nc)*logvar );
    if ((logvar<0) || (AlphaSres>0.7)) {
        AlphaSres = 0.7;
        }
    return Nc/M_PI*AlphaSres;
}

double ComputeSigmaR::alpha_bar_running_pd( double rsq ) { // alphabar = Nc/M_PI * alphas
    double scalefactor = 4.0*alpha_scaling_C2_;
    const double alphas_mu0=2.5;    // mu0/lqcd
    const double alphas_freeze_c=0.2;
    double b0 = (11.0*config::NC - 2.0*config::NF)/3.0;

    double AlphaSres = 4.0*M_PI / ( b0 * std::log(
    std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(scalefactor/(rsq*lambdaqcd*lambdaqcd), 1.0/alphas_freeze_c), alphas_freeze_c)
    ) );
    return Nc/M_PI*AlphaSres;
}

double ComputeSigmaR::alpha_bar_fixed( double rsq ) { // alphabar = Nc/M_PI * alphas
    return 0.190986; // alphas_fixed=0.2
}


// QG term specific running couplings
struct Alphasdata{
    double x01sq,x02sq,x21sq;
};

double ComputeSigmaR::alpha_bar_QG_fixed( void *userdata ) { // alphabar = Nc/M_PI * alphas
    return 0.190986; // alphas_fixed=0.2
}

double ComputeSigmaR::alpha_bar_QG_running_pd( void *userdata ) { // alphabar = Nc/M_PI * alphas
    Alphasdata *dataptr = (Alphasdata*)userdata;
    double rsq=dataptr->x01sq;
    double scalefactor = 4.0*alpha_scaling_C2_;
    const double alphas_mu0=2.5;    // mu0/lqcd
    const double alphas_freeze_c=0.2;
    double b0 = (11.0*config::NC - 2.0*config::NF)/3.0;

    double AlphaSres = 4.0*M_PI / ( b0 * std::log(
    std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(scalefactor/(rsq*lambdaqcd*lambdaqcd), 1.0/alphas_freeze_c), alphas_freeze_c)
    ) );
  return Nc/M_PI*AlphaSres;
}

double ComputeSigmaR::alpha_bar_QG_running_guillaume( void *userdata ) { // alphabar = Nc/M_PI * alphas
    Alphasdata *dataptr = (Alphasdata*)userdata;
    double x01sq=dataptr->x01sq;
    double x02sq=dataptr->x02sq;
    double x21sq=dataptr->x21sq;
    double r_eff_sqr = x01sq * std::pow( x21sq / (x02sq), (x02sq-x21sq)/(x01sq) ); // r = x01 , X = x02 , Y = x21 , Q_{123} = 4C^2 / r_eff_sqr.

    double scalefactor = 4.0*alpha_scaling_C2_;
    const double alphas_mu0=2.5;    // mu0/lqcd
    const double alphas_freeze_c=0.2;
    double b0 = (11.0*config::NC - 2.0*config::NF)/3.0;

    double AlphaSres = 4.0*M_PI / ( b0 * std::log(
    std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(scalefactor/(r_eff_sqr*lambdaqcd*lambdaqcd), 1.0/alphas_freeze_c), alphas_freeze_c)
    ) );
  return Nc/M_PI*AlphaSres;
}

///===========================================================================================
// HIGHER ORDER BK SIGMA3 CORRECTION TERM KERNELS

// resum bk
double ComputeSigmaR::K_resum (double rsq, double x20sq, double x21sq) {
    // returns K_DLA * K_STL
    // double log resummations r = x01 , X = x02 , Y = x21
    double xsq =  std::log(x20sq/rsq) * std::log(x21sq/rsq) ;
    double K_DLA;
    double resummation_alphabar = Alphabar(rsq);
    if (xsq >=0)
        {            
            // argument to the Bessel function is 2sqrt(bar as * xsq)
            double as_x = std::sqrt( resummation_alphabar * xsq );
            K_DLA = gsl_sf_bessel_J1(2.0*as_x) / as_x;
        }
    else // L_xzr L_yzr < 0
        {
            xsq = std::abs(xsq);
            double as_x = std::sqrt( resummation_alphabar * xsq);
            gsl_sf_result res;
            int status = gsl_sf_bessel_I1_e(2.0*as_x, &res);
            
            if (status != GSL_SUCCESS)
            {
                if (std::isnan(x20sq) or std::isnan(x21sq))
                    return 0;	// 0/0, probably z=x or z=y
                cerr << "GSL error " << status <<", result " << res.val << ", as_x=" << as_x << ", xsq=" << xsq << LINEINFO << endl;
                return 0;
            }
            K_DLA = res.val / as_x; 
        }
    if (std::isnan(K_DLA))
        {
            K_DLA =  1;
        }
    // single log resummations
    double K_STL;
    double minxysq = std::min(x20sq,x21sq);
    if (minxysq < 1e-40) minxysq = 1e-40;

    const double A1 = 11.0/12.0;
    K_STL = std::exp( - resummation_alphabar * A1 * std::abs( std::log( config::KSUB * rsq/minxysq ) ) );

    double K_RES = K_DLA*K_STL; // 1601.06598
    return K_RES;
}

double ComputeSigmaR::K_lobk (double rsq, double x20sq, double x21sq) {
    double K_loBK = rsq/(x20sq*x21sq); // 1601.06598
    return K_loBK;
}

///===========================================================================================
// NLO sub scheme qg subtraction terms for different evolution equations.
double ComputeSigmaR::ILNLOqg_subterm_lobk_z2tozero(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit of the full NLO impact factors numerically
    double subterm;
    subterm = ILNLOqg(Q,x,z1, 0,x01sq,x02sq,x21sq);
    return subterm;
}

double ComputeSigmaR::ILNLOqg_subterm_lobk_explicit(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit of the full NLO impact factors explicitly by hand
    double subterm;
    double impact_factor_lo = ILLO(Q,z1,x01sq);
    double P_at_zero = 2;
    double nlo_if_kernel = P_at_zero * 0.5 * (x01sq + x21sq - x02sq) / (x02sq * x21sq); // this is the symmetric half of the LOBK kernel
    double x01 = sqrt(x01sq);
    double x02 = sqrt(x02sq);
    double x21 = sqrt(x21sq);
    double nlo_if_bk_evol_dipoles = Sr(x01, x) - SrTripole(x01, x02, x21, x);
    subterm = impact_factor_lo * nlo_if_kernel * nlo_if_bk_evol_dipoles;
    return subterm;
}

double ComputeSigmaR::ILNLOqg_subterm_resumbk(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit LOBK kernel complemented with the RESUM kernel
    double subterm;
    double K_res = K_resum(x01sq, x02sq, x21sq);
    double lobk_subterm = ILNLOqg_subterm_lobk_explicit(Q,x,z1,z2,x01sq,x02sq,x21sq);
    subterm = K_res * lobk_subterm;
    return subterm;
}

double ComputeSigmaR::ILNLOqg_subterm_kcbk_beuf(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit LOBK kernel complemented with the KCBK kernel from 1401.0313
    double delta012 = std::max(0.0, std::log( std::min(x02sq, x21sq) / (x01sq) ) ); // (166)
    // double shifted_rapidity = helper->rapidity - delta012;
    double scaled_xbj = x * std::exp(delta012); // rapidity shift Y' = Y - Delta012 leads to a scaling of xbj
    if (scaled_xbj > icX0) // Heaviside step with Y_shift > 0 leads to xbj_scaled < 0 and so the cutoff inequality flips.
        return 0;   // Step function in (165)
    
    // Dipoles at shifted rapidity / scaled xbj
    double s02_shifted = Sr(sqrt(x02sq), scaled_xbj);
    double s12_shifted = Sr(sqrt(x21sq), scaled_xbj);
    double s01 = Sr(sqrt(x01sq), x);
    
    double subterm;
    double impact_factor_lo = ILLO(Q,z1,x01sq);
    double P_at_zero = 2;
    double nlo_if_kernel = P_at_zero * 0.5 * (x01sq + x21sq - x02sq) / (x02sq * x21sq); // this is the symmetric half of the LOBK kernel
    subterm = impact_factor_lo * nlo_if_kernel * ( -s02_shifted * s12_shifted + s01);
    return subterm;
}

double ComputeSigmaR::ITNLOqg_subterm_lobk_z2tozero(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit of the full NLO impact factors numerically
    double subterm;
    subterm = ITNLOqg(Q,x,z1,0,x01sq,x02sq,x21sq);
    return subterm;
}

double ComputeSigmaR::ITNLOqg_subterm_lobk_explicit(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit of the full NLO impact factors explicitly by hand
    double subterm;
    double impact_factor_lo = ITLO(Q,z1,x01sq);
    double P_at_zero = 2;
    double nlo_if_kernel = P_at_zero * 0.5 * (x01sq + x21sq - x02sq) / (x02sq * x21sq);
    double x01 = sqrt(x01sq);
    double x02 = sqrt(x02sq);
    double x21 = sqrt(x21sq);
    double nlo_if_bk_evol_dipoles = Sr(x01, x) - SrTripole(x01, x02, x21, x);
    subterm = impact_factor_lo * nlo_if_kernel * nlo_if_bk_evol_dipoles;
    return subterm;
}

double ComputeSigmaR::ITNLOqg_subterm_resumbk(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit LOBK kernel complemented with the RESUM kernel
    double subterm;
    double K_res = K_resum(x01sq, x02sq, x21sq);
    double lobk_subterm = ITNLOqg_subterm_lobk_explicit(Q,x,z1,z2,x01sq,x02sq,x21sq);
    subterm = K_res * lobk_subterm;
    return subterm;
}

double ComputeSigmaR::ITNLOqg_subterm_kcbk_beuf(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit LOBK kernel complemented with the KCBK kernel from 1401.0313
    double delta012 = std::max(0.0, std::log( std::min(x02sq, x21sq) / (x01sq) ) ); // (166)
    // double shifted_rapidity = helper->rapidity - delta012;
    double scaled_xbj = x * std::exp(delta012); // rapidity shift Y' = Y - Delta012 leads to a scaling of xbj
    if (scaled_xbj > icX0) // Heaviside step with Y_shift > 0 leads to xbj_scaled < 0 and so the cutoff inequality flips.
        return 0;   // Step function in (165)
    
    // Dipoles at shifted rapidity / scaled xbj
    double s02_shifted = Sr(sqrt(x02sq), scaled_xbj);
    double s12_shifted = Sr(sqrt(x21sq), scaled_xbj);
    double s01 = Sr(sqrt(x01sq), x);
    
    double subterm;
    double impact_factor_lo = ITLO(Q,z1,x01sq);
    double P_at_zero = 2;
    double nlo_if_kernel = P_at_zero * 0.5 * (x01sq + x21sq - x02sq) / (x02sq * x21sq); // this is the symmetric half of the LOBK kernel
    subterm = impact_factor_lo * nlo_if_kernel * ( -s02_shifted * s12_shifted + s01);
    return subterm;
}


///===========================================================================================
// Computation dipole & data passing helpers
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
///===========================================================================================
///===========================================================================================
// INTEGRATORS for cross sections

/*
// --- L L L --- LO -------- L L L --- LO -------- L L L --- LO -----------
*/
double ComputeSigmaR::ILLO(double Q, double z1, double x01sq) {
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

// MASSLESS
int integrand_ILLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=x01*x01;
    double Xrpdty_lo = Optr->Xrpdty_LO(xbj, Sq(Q));

    double res=(1.0-(Optr->Sr(x01,Xrpdty_lo)))*(Optr->ILLO(Q,z1,x01sq))*x01;
        *f=res;
    return 0;
}

double ComputeSigmaR::LLOp(double Q, double x) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.ComputerPtr=this;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Cuba(cubamethod,ndim,integrand_ILLOp,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

// MASS
int integrand_ILLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double qmass=dataptr->qMass_light;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
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
    Cuba(cubamethod,ndim,integrand_ILLOpMass,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}


/*
// --- T T T --- LO -------- T T T --- LO -------- T T T --- LO -----------
*/

double ComputeSigmaR::ITLO(double Q, double z1, double x01sq) {
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

// MASSLESS
int integrand_ITLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=x01*x01;
    double Xrpdty_lo = Optr->Xrpdty_LO(xbj, Sq(Q));

    double res=(1.0-(Optr->Sr(x01,Xrpdty_lo)))*(Optr->ITLO(Q,z1,x01sq))*x01;
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
    Cuba(cubamethod,ndim,integrand_ITLOp,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

// MASS
int integrand_ITLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double qmass=dataptr->qMass_light;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
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

    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*ef;
    Cuba(cubamethod,ndim,integrand_ITLOpMass,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}



/*
// --- L L L --- NLO -------- L L L --- NLO -------- L L L --- NLO -----------
*/

///*
// NLO
double ComputeSigmaR::ILdip(double Q, double z1, double x01sq) { //ILbeufDIP
    double fac1 = Sq(z1)*Sq(1.0 - z1);
    double facNLO = 4.0*Sq(Q*gsl_sf_bessel_K0(Q*sqrt(z1*(1-z1)*x01sq)));
    return fac1*facNLO;
}

int integrand_ILdip(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) { //integrand_ILbeufDIP
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=Sq(x01);

    double alphabar=Optr->Alphabar(x01sq); //2;
    double alphfac=alphabar*CF/Nc;
    double Xrpdty_lo = Optr->Xrpdty_DIP(xbj, Sq(Q));
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdty_lo);
    double regconst = 5.0/2.0 - Sq(M_PI)/6.0;
    double res;

    res = SKernel*(Optr->ILdip(Q,z1,x01sq))*x01*( alphfac*(0.5*Sq(log(z1/(1-z1))) + regconst ) );
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILdip_z2(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    double res;
    if (z2min > 0.5){ // Check that z2min is not too large since z1 limits are [z2min , 1-z2min]. IF it is too large, return res=0.
        res = 0;
        *f=res;
        return 0;
    }
    double z1=((1.0-z2min)-z2min)*x[0]+z2min;
    double z2=((z1)-z2min)*x[1]+z2min; // z2 integration variable for k1-term AND for the mirror symmetry combined integral!
    double x01=nlodis_config::MAXR*x[2];
    double x01sq=Sq(x01);
    double jac=((1.0-z2min)-z2min)*((z1)-z2min); // jacobiaaniin skaalaus molemmista z-muuttujista!
    double Xrpdt= z2min * X0/z2; // consistent with the qg-terms

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdt); //ok?

    double fivefourths = 5.0/4.0;
    double loopcontribution = 1.0/z2 * ( log(1.0 + z2 / (1.0-z1) ) + log(1.0 - z2 / z1) ) + 1.0/z1 * fivefourths ; // a.k.a k1- or z1-term
    res = 2 * jac * SKernel*(Optr->ILdip(Q,z1,x01sq))*x01*( alphfac*( loopcontribution ) ); // factor of 2 is from the combination of the mirror symmteric q/qbar terms.

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::Bessel0Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    double x01=sqrt(x01sq);
    double x02=sqrt(x02sq);
    double x21=sqrt(x21sq);

    double X3sq = z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq;

    double bessel_innerfun = Q*sqrt(X3sq);
    double facNLO = 0;
    if (bessel_innerfun < 1e-7){
        // cout << "bessel_innerfun = " << bessel_innerfun << " Q " << Q << " X3sq " << X3sq << endl;
        return 0;
    }else{
        facNLO = 4.0*Sq(Q*gsl_sf_bessel_K0( bessel_innerfun ));
    }
    return facNLO*(1-SrTripole(x01,x02,x21,x));
}

double ComputeSigmaR::ILNLOqg(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) { // old ILNLObeufQG
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double fac1 = Sq(z1)*Sq(1.0 - z1);
    double xi 	= z2/(1.0-z1);
    double fun1 = 1+Sq(1-xi);
    double fun2 = x20x21/(x02sq*x21sq);
    double fac2 = 1/x02sq - fun2;
    double fac3 = Sq(xi)*fun2;

    double facNLO1 = Bessel0Tripole(Q, x, z1, z2, x01sq, x02sq, x21sq);
    double facNLO2 = Bessel0Tripole(Q, x, z1, z2, x01sq, 0    , x01sq);

    double res = fac1*((fun1*fac2)*(facNLO1 - facNLO2) + fac3*facNLO1 );
    if(gsl_finite(res)==1){
        return res;
    }else{
        res=0;
    }
    return res;
}

double ComputeSigmaR::ILNLOsigma3(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    double impactfac_lo = ILLO(Q,z1,x01sq);
    double BKkernel = (K_kernel(x01sq,x02sq,x21sq) - 1.0)*2*K_lobk(x01sq,x02sq,x21sq);
    double S012 = SrTripole(sqrt(x01sq),sqrt(x02sq),sqrt(x21sq),x);
    double Sdipole_bkevol = Sr(sqrt(x01sq),x) - S012;

    double res = BKkernel*Sdipole_bkevol*impactfac_lo;
    return res;
}

double ComputeSigmaR::ILNLOqgRisto(double Q, double x, double z0, double z2, double x01sq, double x02sq, double x21sq){ // old ILNLOristo2QG
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);
    double z3 = 1 - z0 - z2;
    double X3sq = z0*z3*x01sq + z0*z2*x02sq + z2*z3*x21sq;

    double fun2 = x20x21/(x02sq*x21sq);

    // FIRST FORM
    double facterm1 = Sq(z3) * ( 2*z0*(z0+z2) + Sq(z2) )/x02sq;
    double facterm2 = Sq(z0) * ( 2*z3*(z2+z3) + Sq(z2) )/x21sq;
    double facterm3 = 2*( (z0+z2)*z0*Sq(z3) + (z2+z3)*z3*Sq(z0) ) * fun2;

    /* // SECOND FORM
    double facterm1 = Sq(z3) * ( 2*z0*(z0+z2) + Sq(z2) ) * ( 1/(x02sq) - fun2 );
    double facterm2 = Sq(z0) * ( 2*z3*(z2+z3) + Sq(z2) ) * ( 1/(x21sq) - fun2 );
    double facterm1expn = Sq(z3) * ( 2*z0*(z0+z2) + Sq(z2) )*(1/x02sq);
    double facterm2expo = Sq(z0) * ( 2*z3*(z2+z3) + Sq(z2) )*(1/x21sq);
    double facterm3 = Sq(z2)*( Sq(z0) + Sq(z3) ) * fun2;
    */

    // double facBesselTrip = Bessel0Tripole(Q, x, z3, z2, x01sq, x02sq, x21sq, X3sq);
    double facBesselTrip = Bessel0Tripole(Q, x, z3, z2, x01sq, x02sq, x21sq);

    double Qbarn = sqrt(z3*(1-z3))*Q;
    double Qbaro = sqrt(z0*(1-z0))*Q;
    double SKernel = 1.0 - Sr(sqrt(x01sq),x);

    ///*
    double facExpn, facExpo;
    if(x02sq/x01sq>1e-8 && x02sq/x01sq<5e2){
        facExpn = gsl_sf_exp(-x02sq/(gsl_sf_exp(M_EULER)*x01sq));
    }else if (x02sq/x01sq<1e-8){
        facExpn = 1;
    }else{ facExpn = 0; }
    if(x21sq/x01sq>1e-8 && x21sq/x01sq<5e2){
        facExpo = gsl_sf_exp(-x21sq/(gsl_sf_exp(M_EULER)*x01sq));
    }else if (x21sq/x01sq<1e-8){
        facExpn = 1;
    }else{ facExpn = 0; }

    double facBesExpn, facBesExpo;
    if(Qbarn*sqrt(x01sq)>1e-8 && Qbarn*sqrt(x01sq)<5e2){
        facBesExpn = 4.0*Sq(Q*gsl_sf_bessel_K0(Qbarn*sqrt(x01sq)))*facExpn*SKernel;
    }else{ facBesExpn = 0; }
    if(Qbaro*sqrt(x01sq)>1e-8 && Qbaro*sqrt(x01sq)<5e2){
        facBesExpo = 4.0*Sq(Q*gsl_sf_bessel_K0(Qbaro*sqrt(x01sq)))*facExpo*SKernel;
    }else{ facBesExpo = 0; }
    //*/

    //double facBesExpn = 4.0*Sq(Q*gsl_sf_bessel_K0(Qbarn*sqrt(x01sq)))*gsl_sf_exp(-x02sq/(gsl_sf_exp(M_EULER)*x01sq))*SKernel;
    //double facBesExpo = 4.0*Sq(Q*gsl_sf_bessel_K0(Qbaro*sqrt(x01sq)))*gsl_sf_exp(-x21sq/(gsl_sf_exp(M_EULER)*x01sq))*SKernel;

    // FIRST FORM
    double term1 = facterm1*( facBesselTrip - facBesExpn );
    double term2 = facterm2*( facBesselTrip - facBesExpo );
    double term3 = facterm3*facBesselTrip;
    double res = term1 + term2 - term3;

    /* //SECOND FORM
    double term1 = facterm1* facBesselTrip - facterm1expn*facBesExpn;
    double term2 = facterm2* facBesselTrip - facterm2expo*facBesExpo;
    double term3 = facterm3*facBesselTrip;
    double res = term1 + term2 + term3;
    */

    return res;
}

int integrand_ILqgunsub(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) { //integrand_ILbeufQGiancu
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }
    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double res =   jac*alphfac*( Optr->ILNLOqg(Q,Xrpdt,z1,z2,x01sq,x02sq,x21sq) )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILsigma3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }
    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double res = jac*alphfac*( Optr->ILNLOsigma3(Q,Xrpdt,z1,0,x01sq,x02sq,x21sq) )/z2*x01*x02; // K_NLO(z_2 = 0)

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgsub(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min >= 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }
    double z1=x[0];
    double z2=(1.0-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min);
    double Xrpdt= z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double ILNLOqg_z2 = 0;
    if (1-z1-z2 > 0){ // heaviside theta front factor on the full NLO term
        ILNLOqg_z2 = (Optr->ILNLOqg(Q,Xrpdt,z1,z2,x01sq,x02sq,x21sq));
    }
    double ILNLOqg_z2_to_0 = Optr->ILNLOqg_subterm(Q,Xrpdt,z1, 0,x01sq,x02sq,x21sq);

    double res = jac*alphfac*(
                    ILNLOqg_z2
                    -
                    ILNLOqg_z2_to_0
                )/(z2)*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgunsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata){
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double x0lim=xbj/X0;
    double z0=(1.0-x0lim)*x[0];
    double z2=(1.0-z0-x0lim)*x[1]+x0lim;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-x0lim)*(1.0-z0-x0lim);
    double Xrpdt=xbj/z2;

    double alphabar=Optr->Alphabar(x01sq); //0.2;
    double alphfac=alphabar*CF/Nc;

    double res;

    res =   jac*alphfac*( Optr->ILNLOqgRisto(Q,Xrpdt,z0,z2,x01sq,x02sq,x21sq) )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata){
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double x0lim=xbj/X0;
    double z0=x[0];
    double z2=(1.0-x0lim)*x[1]+x0lim;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-x0lim);
    double Xrpdt=xbj/z2;

    double alphabar=Optr->Alphabar(x01sq); //0.2;
    double alphfac=alphabar*CF/Nc;

    double res;

    res =   jac*alphfac*(
                    (Optr->ILNLOqgRisto(Q,Xrpdt,z0,z2,x01sq,x02sq,x21sq))*(Optr->heaviside_theta(1-z0-z2))
                    -Optr->ILNLOqgRisto(Q,Xrpdt,z0, 0,x01sq,x02sq,x21sq)
                    )/z2*x01*x02;

        if(gsl_finite(res)==1){
            *f=res;
        }else{
            *f=0;
        }
        return 0;
}

// integrations
double ComputeSigmaR::LNLOdip(double Q, double x) { // old LNLObeufDIP
    double integral, error, prob;
    const int ndim=2;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILdip,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOdip_z2(double Q, double x) { // old LNLObeufDIP
    double integral, error, prob;
    const int ndim=3;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILdip_z2,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOqgunsub(double Q, double x) { // old LNLObeufQGiancu
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILqgunsub,&userdata,&integral,&error,&prob);
    return 2*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOsigma3(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILsigma3,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOqgsub(double Q, double x) { // no old implementation/version. CXY is similar but distinct, sub vs. xbj-sub.
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILqgsub,&userdata,&integral,&error,&prob);
    return 2*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOqgunsubRisto(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILqgunsubRisto,&userdata,&integral,&error,&prob);
    return (2/2)*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral; // TODO a difference of a factor of 2 w.r.t. BEUF formulation was discovered in initial testing
}

double ComputeSigmaR::LNLOqgsubRisto(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILqgsubRisto,&userdata,&integral,&error,&prob);
    return (2/2)*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral; // TODO a difference of a factor of 2 w.r.t. BEUF formulation was discovered in initial testing
}
//*/



/*
// --- T T T --- NLO -------- T T T --- NLO -------- T T T --- NLO -----------
*/

///*
// NLO
double ComputeSigmaR::ITdip(double Q, double z1, double x01sq) {
    double fac1 = z1*(1-z1)*(Sq(z1)+Sq(1.0 - z1));
    double facNLO = Sq(Q*gsl_sf_bessel_K1(Q*sqrt(z1*(1-z1)*x01sq)));
    return fac1*facNLO;
}

int integrand_ITdip(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) { // integrand_ITbeufDIP
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=Sq(x01);

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double Xrpdty_lo = Optr->Xrpdty_DIP(xbj, Sq(Q));
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdty_lo);
    double regconst = 5.0/2.0 - Sq(M_PI)/6.0;
    double res;

    res = SKernel*( Optr->ITdip(Q,z1,x01sq) )*x01*( alphfac*(0.5*Sq(log(z1/(1-z1))) + regconst ) );

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITdip_z2(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    double res;
    if (z2min > 0.5){ // Check that z2min is not too large since z1 limits are [z2min , 1-z2min]. IF it is too large, return res=0.
        res = 0;
        *f=res;
    return 0;
    }
    double z1=((1.0-z2min)-z2min)*x[0]+z2min;
    double z2=((z1)-z2min)*x[1]+z2min; // z2 integration variable for k1-term AND for the mirror symmetry combined integral!
    double x01=nlodis_config::MAXR*x[2];
    double x01sq=Sq(x01);
    double jac=((1.0-z2min)-z2min)*((z1)-z2min); // jacobiaaniin skaalaus molemmista z-muuttujista!
    double Xrpdt= z2min * X0/z2; // consistent with the qg-terms

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdt); //ok?

    double fivefourths = 5.0/4.0;
    double loopcontribution = 1.0/z2 * ( log(1.0 + z2 / (1.0-z1) ) + log(1.0 - z2 / z1) ) + 1.0/z1 * fivefourths ; // a.k.a k1- or z1-term
    res = 2 * jac * SKernel*(Optr->ITdip(Q,z1,x01sq))*x01*( alphfac*( loopcontribution ) ); // factor of 2 is from the combination of the mirror symmteric q/qbar terms.

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::Bessel1Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    double x01=sqrt(x01sq);
    double x02=sqrt(x02sq);
    double x21=sqrt(x21sq);

    double X3sq = z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq;

    double bessel_innerfun = Q*sqrt(X3sq);
    double facNLO = 0;
    if (bessel_innerfun < 1e-7){
        // cout << "bessel_innerfun = " << bessel_innerfun << " Q " << Q << " X3sq " << X3sq << endl;
        return 0;
    }else{
        facNLO = Sq(Q*gsl_sf_bessel_K1(bessel_innerfun));
    }
    double res = facNLO*(1-SrTripole(x01,x02,x21,x));
    return res;
}

double ComputeSigmaR::ITNLOqg(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) { // old ITNLObeufQG
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);
    double X3sq 	= z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq;

    double fac1 	= z1*(1.0 - z1);
    double xi 		= z2/(1.0-z1);
    double fac11	= Sq(z1)+Sq(1-z1);
    double fac12 	= 1+Sq(1-xi);
    double fun2 	= x20x21/(x02sq*x21sq);
    double fac13 	= 1/x02sq - fun2;
    double fac3 	= Sq(xi)*(	fac11*fun2	+ 2*fac1*(1-xi)*x20x21/(x02sq*X3sq)	-	(1-z1)*(1-xi)*(z1+xi-z1*xi)/X3sq	);

    double facNLO1 = Bessel1Tripole(Q, x, z1, z2, x01sq, x02sq, x21sq);
    double facNLO2 = Bessel1Tripole(Q, x, z1, z2, x01sq, 0    , x01sq);

    double res = fac1*((fac11*fac12*fac13)*(facNLO1 - facNLO2) + fac3*facNLO1);
    return res;
}

double ComputeSigmaR::ITNLOsigma3(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    double impactfac_lo = ITLO(Q,z1,x01sq);
    double BKkernel = (K_kernel(x01sq,x02sq,x21sq) - 1.0)*2*K_lobk(x01sq,x02sq,x21sq);
    double S012 = SrTripole(sqrt(x01sq),sqrt(x02sq),sqrt(x21sq),x);
    double Sdipole_bkevol = Sr(sqrt(x01sq),x) - S012;

    double res = BKkernel*Sdipole_bkevol*impactfac_lo;
    return res;
}

double ComputeSigmaR::ITNLOqgRisto(double Q, double x, double z0, double z2, double x01sq, double x02sq, double x21sq){
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);
    double z3 = 1 - z0 - z2;
    double X3sq = z0*z3*x01sq + z0*z2*x02sq + z2*z3*x21sq;

    double fun2 = x20x21/(x02sq*x21sq);

    double factermF = z3/(z0+z2) * ( 1-2*z3*(1-z3) ) * ( 2*z0*(z0+z2) + Sq(z2) )/x02sq;
    double factermG = z0/(z2+z3) * ( 1-2*z0*(1-z0) ) * ( 2*z3*(z2+z3) + Sq(z2) )/x21sq;
    double factermPI = (-2*z0*z3*(
                        (Sq(1-z0)+Sq(z0)) + (Sq(1-z3)+Sq(z3))
                        )*X3sq*fun2
                        +(2*z0*Sq(z2*z3))/(z0+z2)*x20x21/(x02sq)
                        +(2*z3*Sq(z2*z0))/(z2+z3)*x20x21/(x21sq)
                        +z0*z2*z3*(
                            Sq(z0+z2) + Sq(z2+z3) +2*(Sq(z0) + Sq(z3))
                            +Sq(z0*z3)/(Sq(z0+z2))+Sq(z0*z3)/(Sq(z2+z3))
                            -z2*( (z0+z2)/(z2+z3)+ (z2+z3)/(z0+z2) )
                            -( (1-2*z3*(1-z3))*(2*z0*(z0+z2)+Sq(z2)) )/(Sq(z0+z2))
                            -( (1-2*z0*(1-z0))*(2*z3*(z2+z3)+Sq(z2)) )/(Sq(z2+z3))
                            )
                        )/X3sq;

    double facBesselTrip = Bessel1Tripole(Q, x, z3, z2, x01sq, x02sq, x21sq);

    double Qbarn = sqrt(z3*(1-z3))*Q;
    double Qbaro = sqrt(z0*(1-z0))*Q;
    double SKernel = 1.0 - Sr(sqrt(x01sq),x);

    ///*
    double facExpn, facExpo;
    if(x02sq/x01sq>1e-8 && x02sq/x01sq<5e2){
        facExpn = gsl_sf_exp(-x02sq/(gsl_sf_exp(M_EULER)*x01sq));
    }else if (x02sq/x01sq<1e-8){
        facExpn = 1;
    }else{ facExpn = 0; }
    if(x21sq/x01sq>1e-8 && x21sq/x01sq<5e2){
        facExpo = gsl_sf_exp(-x21sq/(gsl_sf_exp(M_EULER)*x01sq));
    }else if (x21sq/x01sq<1e-8){
        facExpn = 1;
    }else{ facExpn = 0; }

    double facBesExpn, facBesExpo;
    if(Qbarn*sqrt(x01sq)>1e-8 && Qbarn*sqrt(x01sq)<5e2){
        facBesExpn = Sq(Q*gsl_sf_bessel_K1(Qbarn*sqrt(x01sq)))*facExpn*SKernel;
    }else{ facBesExpn = 0; }
    if(Qbaro*sqrt(x01sq)>1e-8 && Qbaro*sqrt(x01sq)<5e2){
        facBesExpo = Sq(Q*gsl_sf_bessel_K1(Qbaro*sqrt(x01sq)))*facExpo*SKernel;
    }else{ facBesExpo = 0; }
    //*/
    //double facBesExpn = Sq(Q*gsl_sf_bessel_K1(Qbarn*sqrt(x01sq)))*gsl_sf_exp(-x02sq/(gsl_sf_exp(M_EULER)*x01sq))*SKernel;
    //double facBesExpo = Sq(Q*gsl_sf_bessel_K1(Qbaro*sqrt(x01sq)))*gsl_sf_exp(-x21sq/(gsl_sf_exp(M_EULER)*x01sq))*SKernel;

    double term1 = factermF*( facBesselTrip - facBesExpn );
    double term2 = factermG*( facBesselTrip - facBesExpo );
    double term3 = factermPI*facBesselTrip;
    double res = term1 + term2 + term3;

    return res;
}

int integrand_ITqgunsub(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) { // old integrand_ITbeufQGian
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }
    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double res = jac*alphfac*( Optr->ITNLOqg(Q,Xrpdt,z1,z2,x01sq,x02sq,x21sq) )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITsigma3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }
    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double res =   jac*alphfac*( Optr->ITNLOsigma3(Q,Xrpdt,z1,0,x01sq,x02sq,x21sq) )/z2*x01*x02; // K_NLO(z_2 = 0) kirjoitettu explisiittisesti

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgsub(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min >= 1.0){ // Check that z2min is not too large. IF it is too large, return res=0.
        *f=0;
        return 0;
    }
    double z1=x[0];
    double z2=(1.0-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min);
    double Xrpdt= z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double ITNLOqg_z2 = 0;
    if (1-z1-z2 > 0){ // heaviside theta front factor on the full NLO term
        ITNLOqg_z2 = (Optr->ITNLOqg(Q,Xrpdt,z1,z2,x01sq,x02sq,x21sq));
    }
    double ITNLOqg_z2_to_0 = Optr->ITNLOqg_subterm(Q,Xrpdt,z1,0,x01sq,x02sq,x21sq);

    double res = jac*alphfac*(
                    ITNLOqg_z2
                    -
                    ITNLOqg_z2_to_0
                )/(z2)*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgunsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata){
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double x0lim=xbj/X0;
    double z0=(1.0-x0lim)*x[0];
    double z2=(1.0-z0-x0lim)*x[1]+x0lim;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-x0lim)*(1.0-z0-x0lim);
    double Xrpdt=xbj/z2;

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;

    double res;
    res = jac*alphfac*( Optr->ITNLOqgRisto(Q,Xrpdt,z0,z2,x01sq,x02sq,x21sq) )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata){
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double x0lim=xbj/X0;
    double z0=x[0];
    double z2=(1.0-x0lim)*x[1]+x0lim;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-x0lim);
    double Xrpdt=xbj/z2;

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;

    double res;
    res =   jac*alphfac*(
                    (Optr->ITNLOqgRisto(Q,Xrpdt,z0,z2,x01sq,x02sq,x21sq))*(Optr->heaviside_theta(1-z0-z2))
                    -Optr->ITNLOqgRisto(Q,Xrpdt,z0,0,x01sq,x02sq,x21sq)
                )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}


// integrations
double ComputeSigmaR::TNLOdip(double Q, double x) { // old name: TNLObeufDIP
    double integral, error, prob;
    const int ndim=2;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITdip,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOdip_z2(double Q, double x) { // old name: TNLObeufDIP
    double integral, error, prob;
    const int ndim=3;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITdip_z2,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOqgunsub(double Q, double x) { // old TNLObeufQGiancu
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITqgunsub,&userdata,&integral,&error,&prob);
    return 2*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOsigma3(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITsigma3,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOqgsub(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITqgsub,&userdata,&integral,&error,&prob);
    return 2*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOqgunsubRisto(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITqgunsubRisto,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOqgsubRisto(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITqgsubRisto,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

//*/
