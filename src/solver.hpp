/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2014
 */

#ifndef _NLOBK_SOLVER_H
#define _NLOBK_SOLVER_H

#include "dipole.hpp"
#include <tools/interpolation.hpp>
#include <string>

/* General solver class for the BK equation
 */

class BKSolver
{
    public:
        BKSolver(Dipole* d);    // Constructor takes the dipole amplitude class
        BKSolver();             // Empty constructor, s.t. one can e.g. evaluate alphas
        int Solve(double maxy);	// Solve up to maxy

        double Kernel_lo(double r, double v, double theta);

        double Kernel_nlo(double r, double X, double Y, double X2, double Y2, double z_m_z2);
        double Kernel_nlo_fermion(double r, double X, double Y, double X2, double Y2, double z_m_z2);


        double RapidityDerivative_lo(double r, Interpolator* dipole_interp, double rapidity=-1);
        double RapidityDerivative_nlo(double r, Interpolator* dipole_interp, Interpolator* dipole_interp_s);

        Dipole* GetDipole();
        
        double Alphas(double r);
    
        void SetAlphasScaling(double C2) { alphas_scaling = C2; }

        void SetTmpOutput(std::string fname);
    
        double GetX0() { return x0; }
        void SetX0(double x_) { x0 = x_; }
    
        void SetICX0_nlo_impfac(double x) { icx0_nlo_impfac = x; }
        double GetICX0_nlo_impfac() { return icx0_nlo_impfac; }
    
        void SetICTypicalPartonVirtualityQ0sqr(double x) { icTypicalPartonVirtualityQ0sqr = x; }
        double GetICTypicalPartonVirtualityQ0sqr() { return icTypicalPartonVirtualityQ0sqr; }

	double GetEta0() const { return eta_0; }
	void SetEta0(double e0) { eta_0 = e0; }
    private:
        double alphas_scaling;
        Dipole* dipole;
        std::string tmp_output;         // File which is updated along with the evolution, if empty no temporary results are saved
	double eta_0;	// New eta-BK, eta < eta_0 gives initial condition
    double x0;  // Initial condition refers to xbj, usually=0.01
    double icx0_nlo_impfac; // x0 in the energy conservation requirement, usually =1
    double icTypicalPartonVirtualityQ0sqr;
};






#endif
