/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2015  
 */

#include "nlobk_config.hpp"
#include <sstream>
#include <string>

using namespace std;

using namespace config;

// Default configs
namespace config
{
     bool VERBOSE = false;

     double NC=3;
     double NF=3;

     double LAMBDAQCD = 0.241;
     double LAMBDAQCD2 = LAMBDAQCD*LAMBDAQCD;

     int RINTPOINTS=85;
     int THETAINTPOINTS = 85;
     double INTACCURACY=0.005;
     double MCINTACCURACY = 0.2;
     double MAXR = 10;          // Quite small, only for testing
     double MINR=1e-6;
     unsigned int RPOINTS = 170;

     size_t MCINTPOINTS = 1e7;


     Equation EQUATION = QCD;  

     double DE_SOLVER_STEP = 0.2; // 0.05 paperissa


     double FIXED_AS = 0.2;


     RunningCouplingLO RC_LO = BALITSKY_LO;
     RunningCouplingNLO RC_NLO = PARENT_NLO;

     bool DOUBLELOG_LO_KERNEL = true; // include double log term from the LO kernel
     bool ONLY_DOUBLELOG = false;

     INTEGRATION_METHOD INTMETHOD_NLO = MISER;

     bool LO_BK = false;    // solve only LO BK

     bool FORCE_POSITIVE_N = true;

     bool DNDY=false;

     bool ONLY_NLO = false;

     bool ONLY_LNR = false;
     bool NO_LNR = false;

     bool RESUM_DLOG = false;
     bool RESUM_SINGLE_LOG = false;

     bool NO_K2 = false;

     bool ONLY_RESUM_DLOG = false;
     
     bool ONLY_SUBTRACTION = false;
     
     double KSUB = 1.0;
     
     SINGLELOG_RESUM_RC RESUM_RC = RESUM_RC_PARENT;

     bool ONLY_K1FIN = false;
    
     bool KINEMATICAL_CONSTRAINT = true;
    
    bool EULER_METHOD = true;
}


std::string NLOBK_CONFIG_STRING()
{
    std::stringstream ss;
    
    
    ss << "MC integration method: ";
    if (INTMETHOD_NLO == MISER)
    ss <<"MonteCarlo Miser, points=" << MCINTPOINTS;
    else if (INTMETHOD_NLO == VEGAS)
    ss <<"MonteCarlo Vegas, points=" << MCINTPOINTS;
    else if (INTMETHOD_NLO == MULTIPLE)
    ss << "Multiple integrals (no montecarlo)";
    else
    ss <<"UNKNOWN!";
    ss << ". K1 integration accuracy " << INTACCURACY ;
    ss<< ". LO Kernel RC: ";
    if (RC_LO == FIXED_LO or EQUATION==CONFORMAL_N4)
    ss << " fixed as=" << FIXED_AS;
    else if (RC_LO == SMALLEST_LO)
    ss << " smallest dipole";
    else if (RC_LO == BALITSKY_LO)
    ss << " Balitsky";
    else if (RC_LO == PARENT_LO)
    ss << " Parent dipole";
    else if (RC_LO == PARENT_BETA_LO)
    ss << " Parent dipole, explicit beta";
    else
    ss << " NO STRING IMPLEMENTED!";
    
    ss<< ". NLO Kernel RC: ";
    if (RC_NLO == FIXED_NLO or EQUATION==CONFORMAL_N4)
    ss << " fixed as=" << FIXED_AS;
    else if (RC_NLO == SMALLEST_NLO)
    ss << " smallest dipole";
    else if (RC_NLO  == PARENT_NLO)
    ss << " Parent dipole";
    else
    ss << " NO STRING IMPLEMENTED!";
    
    ss <<". Nc=" << NC << ", Nf=" << NF;
    
    if (EQUATION == QCD)
    {
        if (DOUBLELOG_LO_KERNEL) ss << ". QCD, Double log term in LO kernel included";
        else ss << ". QCD, Double log term in LO kernel NOT included";
    }
    else if (EQUATION == CONFORMAL_QCD) ss << ". Solving for CONFORMAL dipole";
    else if (EQUATION == CONFORMAL_N4) ss << ". Solving in N=4 for CONFORMAL dipole";
    else ss << ". UNKNOWN EQUATION!!";
    
    
    if (FORCE_POSITIVE_N)
    ss << ". Amplitude is limited to [0,1].";
    else
    ss << ". Amplitude is not limited!";
    
    ss << endl;
    //BKSolver sol;
    //ss << "# Alphas(r=1 GeV^-1) = " << sol.Alphas(1) << endl;
    ss << "# Order: ";
    if (LO_BK)
    ss <<"LO";
    else
    ss << "NLO";
    if (config::ONLY_NLO) ss << ", keeping only NLO terms";
    if (config::RESUM_DLOG)
    {
        ss << endl;
        ss << "# Resumming double log";
    }
    if (config::RESUM_SINGLE_LOG)
    {
        ss << endl;
        ss << "# Resumming single log, K_sub=" << config::KSUB;
        if (config::RESUM_RC == RESUM_RC_PARENT) ss << " resum rc: parent";
        else if (config::RESUM_RC == RESUM_RC_SMALLEST) ss << " resum rc: smallest";
        else if (config::RESUM_RC == RESUM_RC_BALITSKY) ss << " resum rc: balitsky";
        ss << endl;
    }
    
    if (config::ONLY_SUBTRACTION)
    ss << endl << "# Only including the subtraction term" << endl;
    
    if (config::ONLY_K1FIN)
    ss << endl << "# Only including K1fin part of K1" << endl;
    
    if (config::NO_K2)
    {
        ss << endl << "# Not including K2 and Kf" << endl;
    }
    
    if (config::KINEMATICAL_CONSTRAINT)
        ss << endl << "# Kinematical constraint included" << endl;
    return ss.str();
}

