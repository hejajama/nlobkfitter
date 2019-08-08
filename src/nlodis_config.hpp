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

    enum SubtractionScheme
    {
        SUBTRACTED,     // \sigma_NLO = \sigma_LO + \sigma_qg,sub
        UNSUBTRACTED    // \sigma_NLO = \sigma_IC + \sigma_qg,unsub
    };
    extern SubtractionScheme SUB_SCHEME;

    enum SubSchemeTermKernel
    {
        SUBTERM_LOBK_Z2TOZERO,
        SUBTERM_LOBK_EXPLICIT,
        SUBTERM_RESUM,
        SUBTERM_KCBK_BEUF,
        SUBTERM_TRBK_EDMOND
    };
    extern SubSchemeTermKernel SUB_TERM_KERNEL;

    enum GluonZ2Minimum
    {
        Z2SIMPLE,   // z2min = xbj/x0
        Z2IMPROVED  // z2min = (xbj/x0)*(Q0^2/Q^2)
    };
    extern GluonZ2Minimum Z2MINIMUM;

    enum RunningCouplingDIS
    {
        DIS_RC_FIXED,
        DIS_RC_PARENT,      // Parent dipole where all beta terms are included
		DIS_RC_GUILLAUME    // 1708.06557
    };
    extern RunningCouplingDIS RC_DIS;

    enum TargetRapidityBKRhoPresc
    {
        TRBK_RHO_X_R,
        TRBK_RHO_MAX_X_Y_R
    };
    extern TargetRapidityBKRhoPresc TRBK_RHO_PRESC;
}

#endif