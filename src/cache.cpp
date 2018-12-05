
/*
 * Cache BK solution s.t. the code does not need to solve it again
 * if only sigma_r changes
 *
 * The idea is that the previous bk solution is always saved on the disk
 *
 * NOTE: THIS IS VERY DANGEROUS PIECE OF CODE
 * If MNINUIT is not complied in a single-thread mode, crazy things happen!!!
 */

#include "cache.hpp"
#include <vector>
#include <Minuit2/MnUserParameterState.h>
#include <cstdlib>
#include <cmath>

BKCache bkcache;

using namespace ROOT::Minuit2;

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
