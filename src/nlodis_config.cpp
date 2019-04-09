/*
 * LCPT NLO DIS fitter
 * Henri Hänninen <henri.j.hanninen@student.jyu.fi>, 2017-
 */

#include "nlodis_config.hpp"


// Default configs
namespace nlodis_config
{
    int CUBA_MAXEVAL;
    double CUBA_EPSREL;
    bool VERBOSE = false;
    bool PRINTDATA = false;

    double MAXR;
    double MINR;

    RunningCouplingDIS RC_DIS;
}