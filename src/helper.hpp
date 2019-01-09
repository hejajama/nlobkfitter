/*
 * Some helper stuff for BK DIS fits, e.g. find optimal sigma02
 */

#ifndef _CACHE_H
#define _CACHE_H

#include <Minuit2/MnUserParameterState.h>
#include <vector>
// Optimize sigma02 (the coefficient which multiplies SigmaComputer.SigmarLOmass
// Returns optimal sigma02 and total chi^2
std::vector<double> FindOptimalSigma02(std::vector<double> expdata, std::vector<double> experr,
                         std::vector<double> thdata);



#endif
