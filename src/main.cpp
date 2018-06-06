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
{
    errors_mmyiss++;
    std::cerr << file << ":"<< line <<": Error " << errors_mmyiss << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}

int main()
{
    gsl_set_error_handler(&ErrHandlerCustom);
    //gsl_set_error_handler_off();

        config::RC_LO = config::BALITSKY_LO; //config::GUILLAUME_LO; // Balitsky running coupling for LO kernel
        config::RESUM_RC = config::RESUM_RC_PARENT; // Parent dipole in the resummation
        config::RESUM_DLOG = true; // Resum doulbe logs
        config::RESUM_SINGLE_LOG = true; // Resum single logs
        config::LO_BK = false;  // Solve LO BK with running coupling, overrides RESUM settings
        config::KSUB = 0.65;  // Optimal value for K_sub
        config::NO_K2 = true;  // Do not include numerically demanding full NLO part
        config::INTACCURACY = 0.02;//0.02;
        config::MINR = 1e-5;
        config::MAXR = 25;
        config::RPOINTS = 100;
        config::DE_SOLVER_STEP = 0.2; // Euler method probably requires smaller step!
		config::DNDY=false;
        //sigmar_config::maxy = 5.2;

        // If want to use kinematical constraint in the LO equation
        config::EULER_METHOD = false;        // Kinematical constraint requires this
        config::KINEMATICAL_CONSTRAINT = false;

        // Constants
        config::NF=3;   // Only light quarks
        config::LAMBDAQCD = 0.241;


    Data data;
    data.SetMinQsqr(1.0);
    data.SetMaxQsqr(50);
    data.SetMaxX(0.01); // NLO, ic at xbj=1e-1;  TODO Is this necessary?

    // Add datafiles, if 2nd parameter=CHARM, then this is only charmdata
    data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);
    //data.LoadData("data/hera_combined_sigmar_cc.txt", CHARM); // charm data


    MnUserParameters parameters;
	 //parameters.SetPrecision(0.001); 
      // Constants
        //parameters.Add("anomalous_dimension", 1.0);

      // Fit parameters, second value is starting value, second is uncertainty
        /*
        // MASSLESS -- LO LO LO
          parameters.Add("qs0sqr", 0.0564335 , 0.04 );
          parameters.Add("fitsigma0", 42.0015 , 5.0); // 1mb = 2.568 GeV² // (2.568)*16.36
          parameters.Add("alphascalingC2", 7.2 , 1.0 );
          parameters.Add("e_c", 18.9692 , 0.5 );
        */

        ///*
        // MASSLESS -- NLO NLO NLO

          // MV
    	  
          /*parameters.Add("qs0sqr", 0.104, 0.1 );
          parameters.Add("fitsigma0", 48.304, 0.1 ); // 1mb = 2.568 GeV² // (2.568)*16.36
          parameters.Add("alphascalingC2", 14.5 ,1);
          parameters.Add("e_c", 1.0 );
          parameters.Add("anomalous_dimension", 1.0 );
		  parameters.Add("initialconditionX0", 0.01 );
         */
		  // MV for resummed
		  parameters.Add("qs0sqr", 0.1833036828274, 0.1);
		parameters.Add("fitsigma0", 26.24896319135, 2 ); // 1mb = 2.568 GeV² // (2.568)*16.36
          parameters.Add("alphascalingC2", 0.05862105564887, 0.03);
          parameters.Add("e_c", 1.5884630748, 0.5);
          parameters.Add("anomalous_dimension", 1.0 );
		  parameters.Add("initialconditionX0", 0.01 );

          /*
          parameters.Add("qs0sqr", 0.0477335 , 0.04 );
          parameters.Add("fitsigma0", 91.5015 , 5.0); // 1mb = 2.568 GeV² // (2.568)*16.36
          parameters.Add("alphascalingC2", 10.39 , 3.0 ); //alphas_scaling = exp(-2.0*0.5772);
          parameters.Add("e_c", 1.0 );
          parameters.Add("anomalous_dimension", 1.0 );
          parameters.Add("initialconditionX0", 0.01 );
          */

          
          // MVe
		  /*
          parameters.Add("qs0sqr", 0.0564335 , 0.04 );
          parameters.Add("e_c", 18.9692 , 0.5 );
          parameters.Add("fitsigma0", 42.0015 , 5.0); // 1mb = 2.568 GeV² // (2.568)*16.36
          parameters.Add("alphascalingC2", 7.2 , 1.0 );
         
  parameters.Add("anomalous_dimension", 1.0 );
          parameters.Add("initialconditionX0", 0.01 );
*/
        //*/

         /*
        // MASS
          parameters.Add("qs0sqr", 0.06 , 0.04 );
          parameters.Add("e_c", 18.9 , 0.5 );
          parameters.Add("fitsigma0", (2.568)*16.36 , 5.0); // 1mb = 2.568 GeV² // (2.568)*16.36
          parameters.Add("alphascalingC2", 7.2 , 1.0 );
         */

      /*// Set limits
        parameters.SetLowerLimit("A_g", 0);
        parameters.SetUpperLimit("mu_0", 1.43);
      */

    NLODISFitter fitter(parameters);
    fitter.AddDataset(data);
    fitter.SetNLO(false);

    cout << "=== Initial parameters ===" << endl;
    cout << parameters << endl;
    cout << "=== Starting fit ===" << endl;

    // MnMinimize: use MIGRAD, if it fails, fall back to SIMPLEX
    //MnSimplex fit(fitter,parameters);
    MnMinimize fit(fitter, parameters);
    //MnMigrad fit(fitter, parameters);
    //MnScan fit(fitter, parameters);

    // minimize
    FunctionMinimum min = fit();
    // output
    std::cout<<"minimum: "<<min<<std::endl;

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
