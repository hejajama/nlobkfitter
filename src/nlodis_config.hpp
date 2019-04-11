/*
 * LCPT NLO DIS fitter
 * Henri Hänninen <henri.j.hanninen@student.jyu.fi>, 2017-
 */

// Configs

#ifndef _NLODIS_CONFIG_HPP
#define _NLODIS_CONFIG_HPP


namespace nlodis_config
{
    extern int      CUBA_MAXEVAL;
    extern double   CUBA_EPSREL;
    extern bool     VERBOSE;
    extern bool     PRINTDATA;

    extern double MAXR;
    extern double MINR;

    enum RunningCouplingDIS
    {
        DIS_RC_FIXED,
        DIS_RC_PARENT,      // Parent dipole where all beta terms are included
		DIS_RC_GUILLAUME    // 1708.06557
    };
    extern RunningCouplingDIS RC_DIS;

}

#endif