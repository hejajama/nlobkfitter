/*
 * Cache BK solution s.t. the code does not need to solve it again
 * if only sigma_r changes
 *
 * The idea is that the previous bk solution is always saved on the disk
 *
 * NOTE: THIS IS VERY DANGEROUS PIECE OF CODE
 * If MNINUIT is not complied in a single-thread mode, crazy things happen!!!
 */

#ifndef _CACHE_H
#define _CACHE_H

#include <Minuit2/MnUserParameterState.h>
#include <vector>

struct BKCache
{
    std::vector <double> params;     // Parameters that correspond to previous bk solution
};

extern BKCache bkcache;

bool IsResultCached(std::vector<double> values, ROOT::Minuit2::MnUserParameters parameters);
#endif
