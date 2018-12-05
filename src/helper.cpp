
/*
 * Cache BK solution s.t. the code does not need to solve it again
 * if only sigma_r changes
 *
 * The idea is that the previous bk solution is always saved on the disk
 *
 * NOTE: THIS IS VERY DANGEROUS PIECE OF CODE
 * If MNINUIT is not complied in a single-thread mode, crazy things happen!!!
 */

#include "helper.hpp"
#include <vector>
#include <iostream>
#include <Minuit2/MnUserParameterState.h>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

BKCache bkcache;

using namespace std;

using namespace ROOT::Minuit2;



// Optimize sigma02 (the coefficient which multiplies SigmaComputer.SigmarLOmass
// Returns optimal sigma02

struct sigma02helper
{
    std::vector<double> expdata;
    std::vector<double> experr;
    std::vector<double> thdata;
    
};

double minimiser_helper_chisqr(double sigma02, void* p)
{
    sigma02helper *par = (sigma02helper*)p;
    
    double sum=0;
    for (unsigned int i=0; i< par->expdata.size(); i++)
    {
        sum += std::pow((sigma02*par->thdata[i]-par->expdata[i])/par->experr[i], 2.0)/par->expdata.size();
    }
    //cout << sigma02 << " " << sum << endl;
    return sum;
}

std::vector<double> FindOptimalSigma02(std::vector<double> expdata, std::vector<double> experr,
                          std::vector<double> thdata)
{
    sigma02helper par;
    par.expdata = expdata;
    par.experr = experr;
    par.thdata = thdata;
    
    if (expdata.size() != experr.size() or expdata.size() != thdata.size())
    {
        cerr << "Optimizer did not get vectors with same size???" << endl;
        exit(1);
    }
    
    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    double min = 30.0;
    double lower = 1; double upper = 150;;
    gsl_function F;
    
    F.function = &minimiser_helper_chisqr;
    F.params = &par;
    
    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T);
    gsl_min_fminimizer_set (s, &F, min, lower, upper);
    
    
    do
    {
        iter++;
        status = gsl_min_fminimizer_iterate (s);
        
        min = gsl_min_fminimizer_x_minimum (s);
        lower = gsl_min_fminimizer_x_lower (s);
        upper = gsl_min_fminimizer_x_upper (s);
        
        status = gsl_min_test_interval (lower, upper, 0, 0.0001);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_min_fminimizer_free (s);
    
    //cout << "Minimizer ended, iterations " << iter << " minimum " << min << endl;
    
    std::vector<double> result;
    result.push_back(min);
    result.push_back(minimiser_helper_chisqr(min, &par));
    return result;
}



bool IsClose(double a, double b, double relative=0.0000001)
{
	if (a < 0 and b > 0)
		return false;
	if (b<0 and a>0)
		return false;

    if (std::abs(a-b)/std::min( std::abs(a), std::abs(b))<relative)
        return true;
    else
        return false;
}

bool IsResultCached(std::vector<double> pars, MnUserParameters parameters)
{
    // Check if the BK parameters are the same as in the cached solution
    if (pars.size() != bkcache.params.size())
        return false;
    
    double qs0sqr       = pars[ parameters.Index("qs0sqr")];
    double e_c          = pars[ parameters.Index("e_c")];
    double alphas_scaling       = pars[ parameters.Index("alphascalingC2")];
    double anomalous_dimension  = pars[ parameters.Index("anomalous_dimension")];
    double initialconditionX0  = pars[ parameters.Index("initialconditionX0")];
    
    double qs0sqr_cache       = bkcache.params[ parameters.Index("qs0sqr")];
    double e_c_cache          = bkcache.params[ parameters.Index("e_c")];
    double alphas_scaling_cache       = bkcache.params[ parameters.Index("alphascalingC2")];
    double anomalous_dimension_cache  = bkcache.params[ parameters.Index("anomalous_dimension")];
    double initialconditionX0_cache  = bkcache.params[ parameters.Index("initialconditionX0")];
    
    if (IsClose(qs0sqr, qs0sqr_cache) and IsClose(e_c, e_c_cache) and
        IsClose(alphas_scaling, alphas_scaling_cache) and IsClose(anomalous_dimension, anomalous_dimension_cache)
        and IsClose(initialconditionX0, initialconditionX0_cache))
        return true;
    else
        return false;
    
}
