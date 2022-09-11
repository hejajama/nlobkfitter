#include "mv-nucl.hpp"
#include <sstream>
#include <string>
#include <cmath>
#include "nlobk_config.hpp"
#include <iostream>
#include "solver.hpp"
#include <tools/tools.hpp>

using std::cerr; using std::endl;

double MVnuc::DipoleAmplitude(double r, double b)
{
	if (b>1e-10)
		cerr << "Impact parameter dependent BK evolution is not supported!" << LINEINFO << endl;
	if (r < 1e-30)
		return 0;
	const double e = 2.7182818;
	///TODO: some algorithm to determine small r, e.g. when one has to linearize
    if (r < 2e-6)   ///NOTE: factor 1/4 "correctly", not as in AAMS paper
            return A*T_A(impact_b,A)*sigma0/2*std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0
            * std::log( 1.0/(r*lambdaqcd) + ec*e) ;
    return 1.0 - std::exp(-A*T_A(impact_b,A)*sigma0/2*std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0
            * std::log( 1.0/(r*lambdaqcd) + ec*e) );
}


void MVnuc::SetQsqr(double qsqr)
{
	qs0sqr=qsqr;
}

void MVnuc::SetAnomalousDimension(double gamma_)
{
	anomalous_dimension=gamma_;
}

void MVnuc::SetLambdaQcd(double lambda)
{
	lambdaqcd=lambda;
}

void MVnuc::SetE(double ec_)
{
	ec=ec_;
}

void MVnuc::setImpactParb(double b_)
{
	impact_b=b_;
}

void MVnuc::setSigma0(double sigma0_)
{
	sigma0=sigma0_;
}

void MVnuc::setA(int A_)
{
	A=A_;
	InitializeWSDistribution(A);
}

double MVnuc::GetE()
{
    return ec;
}

std::string MVnuc::GetString()
{
	std::stringstream ss;
	ss << "Optical glauber integrated MVnuc model, Q_s0^2 = " << qs0sqr << " GeV^2, \\gamma = " << anomalous_dimension
		<< ", coefficient of E inside Log is " << ec
		<< ", x0=" << x0 << ", \\Lambda_QCD = " << lambdaqcd << " GeV";
	return ss.str();
}

/*
 * Set some reasonable parameters
 * That is, MV1 model for nucleus
 */
MVnuc::MVnuc()
{
	qs0sqr = 0.10;
	x0=0.01;
	ec=1.0;
	impact_b = 0;
	sigma0 = 10;
	// A = 100;
	lambdaqcd=config::LAMBDAQCD;
	anomalous_dimension=1;
}
