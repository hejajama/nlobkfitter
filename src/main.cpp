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

#include "gitsha1.h"

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
	cout << "# NLOBKDISFitter, git version " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << " main build " << __DATE__  << " " << __TIME__     << endl;
    gsl_set_error_handler(&ErrHandlerCustom);
    //gsl_set_error_handler_off();

        config::RC_LO = config::GUILLAUME_LO;// FIXED_,PARENT_,PARENT_BETA_,SMALLEST_,BALITSKY_,FRAC_,GUILLAUME_,
        config::RESUM_RC = config::RESUM_RC_GUILLAUME; // _BALITSKY,_PARENT,_SMALLEST,_GUILLAUME,
        config::RESUM_DLOG = false;; // Resum doulbe logs
        config::RESUM_SINGLE_LOG = false; // Resum single logs
        config::LO_BK = false; // Solve LO BK with running coupling, overrides RESUM settings
        config::KSUB = 0.65;  // Optimal value for K_sub
        config::NO_K2 = true;  // Do not include numerically demanding full NLO part
        config::INTACCURACY = 0.015;//0.02;
        config::MINR = 1e-5;
        config::MAXR = 50;
        config::RPOINTS = 100;
        config::DE_SOLVER_STEP = 0.050;// Euler method probably requires smaller step!
		config::DNDY=false;
        //sigmar_config::maxy = 5.2;

        // If want to use kinematical constraint in the LO equation
        config::EULER_METHOD = true;;    // Kinematical constraint requires this
        config::KINEMATICAL_CONSTRAINT = config::KC_NONE;;

        // Constants
        config::NF=3;   // Only light quarks
        config::LAMBDAQCD = 0.241;


    Data data;
    data.SetMinQsqr(0.75);
    data.SetMaxQsqr(50);
    data.SetMaxX(0.01); // NLO, ic at xbj=1e-1;  TODO Is this necessary?

    // Add datafiles, if 2nd parameter=CHARM, then this is only charmdata
	data.LoadData("./data/hera_combined_sigmar.txt", TOTAL);
//    data.LoadData("./data/light_quark_f2/hera_I_combined_eplus_lightq", TOTAL);
    //data.LoadData("data/hera_combined_sigmar_cc.txt", CHARM); // charm data

	MnMachinePrecision precision;
////	precision.SetPrecision(0.01);
	MnUserParameters parameters;
	 //parameters.SetPrecision(0.001); 
      // Constants
    
    parameters.Add("qs0sqr", 0.0484000000000, 0.4);
    parameters.Add("alphascalingC2", 14.4      , 0.70);
    parameters.Add("e_c", 21.5, 1.0);
    parameters.Add("anomalous_dimension", 1.0 );
	
    // eta_0 starint value for the k^- evolution
	parameters.Add("eta0", std::log(1/0.01));
    
    // Constants
    parameters.Add("initialconditionX0", 1 );
    parameters.Add("initialconditionY0", 0 ); // What rapidity is our IC
    parameters.Add("icTypicalPartonVirtualityQ0sqr", 1.0);
    
      /*// Set limits
        parameters.SetLowerLimit("A_g", 0);
        parameters.SetUpperLimit("mu_0", 1.43);
      */

    NLODISFitter fitter(parameters);
    fitter.AddDataset(data);
    fitter.SetNLO(false);
    
    cout <<"=== BK solver setup: ===" << endl;
    cout << NLOBK_CONFIG_STRING()<<endl;

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
